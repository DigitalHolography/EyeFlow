function [maskArtery, maskVein, scoreMaskArtery, scoreMaskVein] = createMasksSegmentationNet(M0_ff, net, maskArtery)
% createMasksSegmentationNet Create artery and vein masks using a SegmentationNet
%   [maskArtery, maskVein, scoreMaskArtery, scoreMaskVein, R] = createMasksSegmentationNetInternal(M0_ff, M0_ff_img, net, maskArtery);
%   Inputs:
%       M0_ff: 3D array (numX, numY, Nt) of flow-encoded images over time
%       net: trained SegmentationNet model
%       maskArtery: 2D binary array (numX, numY) pre-mask of artery locations
%   Outputs:
%       maskArtery: 2D binary array (numX, numY) artery mask
%       maskVein: 2D binary array (numX, numY) vein mask
%       scoreMaskArtery: 2D array (numX, numY) of artery mask scores
%       scoreMaskVein: 2D array (numX, numY) of vein mask scores

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
mask_params = params.json.Mask;

M0_ff_img = squeeze(mean(M0_ff, 3));
[numX, numY] = size(M0_ff_img);

% diasys_diff = [];
% R = [];
% exportVideos = params.exportVideos;

if ~any(maskArtery, "all")
    error("Given Pre Mask Artery is completely empty.");
end

% if the correlation map is used by the model, compute it
if mask_params.AVCorrelationSegmentationNet
    fprintf("Compute correlation for artery/vein segmentation\n");

    signal = sum(M0_ff .* maskArtery, [1 2], 'omitnan');
    signal = signal ./ nnz(maskArtery);

    outlier_frames_mask = isoutlier(signal, "movmedian", 5, "ThresholdFactor", 2);
    video = interpolateOutlierFrames(M0_ff, outlier_frames_mask);

    signal = sum(video .* maskArtery, [1 2], 'omitnan');
    signal = signal ./ nnz(maskArtery);

    % compute local-to-average signal wave zero-lag correlation
    signal_centered = signal - mean(signal, 3, 'omitnan');
    video_centered = video - mean(M0_ff, 'all', 'omitnan');
    R = mean(video_centered .* signal_centered, 3) ./ (std((video_centered), [], 'all', 'omitnan') * std(signal_centered, [], 3));
    saveMaskImage(R, 'all_15_Correlation.png', isStep = true)
end

% if the systolic and diastolic frames are used by the model, compute them
if mask_params.AVDiasysSegmentationNet
    fprintf("Compute diastolic and stystolic frames for artery/vein segmentation\n");

    [M0_Systole_img, M0_Diastole_img] = compute_diasys(M0_ff, maskArtery, 'mask');
    saveMaskImage(rescale(M0_Systole_img), 'artery_21_systole_img.png', isStep = true)
    saveMaskImage(rescale(M0_Diastole_img), 'vein_21_diastole_img.png', isStep = true)

    diasys_diff = M0_Systole_img - M0_Diastole_img;

    RGBdiasys = labDuoImage(rescale(M0_ff_img), diasys_diff);
    saveMaskImage(RGBdiasys, 'vessel_21_diasys_diff.png', isStep = true);

    M0_Diastole_img = imresize(rescale(M0_Diastole_img), [512, 512]);
    M0_Systole_img = imresize(rescale(M0_Systole_img), [512, 512]);
    diasys_diff = imresize(rescale(diasys_diff), [512, 512]);

end

M0 = imresize(rescale(M0_ff_img), [512, 512]);

if mask_params.AVCorrelationSegmentationNet
    R = imresize(rescale(R), [512, 512]);

    if mask_params.AVDiasysSegmentationNet
        input = cat(3, M0, R, diasys_diff);
    else
        input = cat(3, M0, M0, R);
    end

elseif mask_params.AVDiasysSegmentationNet
    input = cat(3, M0, M0_Diastole_img, M0_Systole_img);
end

ToolBox.Output.add("NetworkInput", input, h5path = "/Maps/NetworkInput");
output = runAVInference(net, input);

[~, argmax] = max(output, [], 3);

onehot = multi2onehot(argmax, 3);

maskArtery = onehot(:, :, 1);
maskVein = onehot(:, :, 2);

scoreMaskArtery = sum(maskArtery .* output(:, :, 2), [1, 2]) / nnz(maskArtery);
scoreMaskVein = sum(maskVein .* output(:, :, 3), [1, 2]) / nnz(maskVein);

maskArtery = imresize(maskArtery, [numX, numY], "nearest");
maskVein = imresize(maskVein, [numX, numY], "nearest");

maskArtery = logical(maskArtery);
maskVein = logical(maskVein);

% if exportVideos && ~isempty(diasys_diff)
%     RGB_diasys_video = labDuoVideo(M0_ff, diasys_diff);
%     writeGifonDisc(mat2gray(RGB_diasys_video), 'diasys.gif');
% elseif exportVideos && ~isempty(R)
%     RGB_corr_video = labDuoVideo(M0_ff, R);
%     writeGifonDisc(mat2gray(RGB_corr_video), 'correlation.gif');
% end

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

function output = runAVInference(model_struct, input)

if model_struct.use_onnx

    if canUseGPU
        device = 'cuda';
    else
        device = 'cpu';
    end

    if isa(model_struct.onnx_model, 'dlnetwork')
        % For dlnetwork objects
        input_dl = dlarray(input, 'SSCB'); % Convert to dlarray
        output_dl = predict(model_struct.onnx_model, input_dl);
        output = extractdata(output_dl);
    elseif isa(model_struct.onnx_model, 'DAGNetwork')
        % For DAGNetwork objects
        output = predict(model_struct.onnx_model, input, 'ExecutionEnvironment', device);
    end

elseif model_struct.use_pytorch
    np = py.importlib.import_module('numpy');
    torch = py.importlib.import_module('torch');

    if torch.cuda.is_available()
        device = 'cuda';
    else
        device = 'cpu';
    end

    % MATLAB HWC → NumPy CHW
    input_np = np.array(permute(input, [3 1 2]), 'float32');

    % To tensor
    t = torch.tensor(input_np).unsqueeze(int32(0)).to(device); % BCHW

    % Forward pass
    model = model_struct.py_model.to(device);
    out = model(t);
    out_np = out.cpu().detach().numpy();

    % NumPy BCHW → MATLAB HWC
    out_mat = squeeze(single(out_np)); % C × H × W
    output = gather(out_mat);

    output = permute(output, [2 3 1]); % CHW to HWC
end

end
