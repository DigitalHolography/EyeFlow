function [maskArtery, maskVein] = createMasksHolonet(M0_ff_video, M0_ff_img, maskVesselness)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
mask_params = params.json.Mask;

% Downloads all huggingface models once
if ~isfolder('Models')
    mkdir('Models');
end

if ~isfile('Models\iternet5_av_diasys.onnx')
    % Download the model from Hugging Face
    url = 'https://huggingface.co/DigitalHolography/iternet5_av_diasys/resolve/main/iternet5_av_diasys';
    websave('Models\iternet5_av_diasys.onnx', url);
end


    function onehot = multi2onehot(x, axis)
    % Converts a label mask x to a binary one-hot encoding with 2 classes
    % axis: dimension along which to stack (default 3)
    
    if nargin < 2
        axis = 3; % Default stacking axis
    end
    
    % Create binary masks
    class1 = (x == 2) | (x == 4);
    class2 = (x == 3) | (x == 4);
    
    % Stack along specified axis
    switch axis
        case 1
            onehot = cat(1, class1, class2);
        case 2
            onehot = cat(2, class1, class2);
        case 3
            onehot = cat(3, class1, class2);
        otherwise
            error('Unsupported axis value: %d', axis);
    end
    
    % Convert to logical or double as needed
    onehot = double(onehot);
    end

[Nx, Ny] = size(M0_ff_img);

% if the correlation map is used by the model, compute it
if mask_params.CorrelationHolonet
    signal = sum(M0_ff_video .* maskVesselness, [1 2], 'omitnan');
    signal = signal ./ nnz(maskVesselness);
    
    outlier_frames_mask = isoutlier(signal, "movmedian", 5, "ThresholdFactor", 2);
    video = interpolateOutlierFrames(M0_ff_video, outlier_frames_mask);
    
    signal = sum(video .* maskVesselness, [1 2], 'omitnan');
    signal = signal ./ nnz(maskVesselness);
    
    % compute local-to-average signal wave zero-lag correlation
    signal_centered = signal - mean(signal, 3, 'omitnan');
    video_centered = video - mean(M0_ff_video, 'all', 'omitnan');
    R = mean(video_centered .* signal_centered, 3) ./ (std((video_centered), [], 'all', 'omitnan') * std(signal_centered, [], 3));
    saveImage(R_VascularSignal, 'all_15_Correlation.png', isStep = true)
end

% if the systolic and diastolic frames are used by the model, compute them
if mask_params.DiaSysHolonet
    [M0_Systole_img, M0_Diastole_img, ~] = compute_diasys(M0_ff_video, maskVesselness, 'mask');
    saveImage(rescale(M0_Systole_img), 'artery_20_systole_img.png', isStep = true)
    saveImage(rescale(M0_Diastole_img), 'vein_20_diastole_img.png', isStep = true)
    M0_Diastole_img = imresize(rescale(M0_Diastole_img), [512, 512]);
    M0_Systole_img = imresize(rescale(M0_Systole_img), [512, 512]);
end

M0 = imresize(rescale(M0_ff_img), [512,512]);

if mask_params.CorrelationHolonet
    R = imresize(rescale(R), [512, 512]);
    if mask_params.DiaSysHolonet
        net = importONNXNetwork('Models\iternet_5_av_sys_corr.onnx');
        
        input = cat(3, M0, M0_Systole_img, M0_Diastole_img, R);
        output = predict(net, input);
    else
        net = importONNXNetwork('Models\iternet_5_av_corr.onnx');
        input = cat(3, M0, M0, R);
        output = predict(net, input);
    end
elseif mask_params.DiaSysHolonet
    net = importONNXNetwork('Models\iternet5_av_diasys.onnx');
    input = cat(3, M0, M0_Diastole_img, M0_Systole_img);
    output = predict(net, input);
end

[~, argmax] = max(output, [], 3);

onehot = multi2onehot(argmax, 3);

maskArtery = onehot(:,:,1);
maskVein = onehot(:,:,2);

maskArtery = imresize(maskArtery, [Nx, Ny], "nearest");
maskVein = imresize(maskVein, [Nx, Ny], "nearest");

end


