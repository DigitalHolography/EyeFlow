function stabilizeVideo(refImgPath, targetVideo, outGifPath)
    refImg = safeLoadGray(refImgPath);
    refImg = imresize(refImg, [512 512]);
    refImg = flat_field_correction(refImg, 10, 0, 'gaussianBlur');

    v = VideoReader(targetVideo);
    nFrames = floor(v.Duration * v.FrameRate) / 5;
    fprintf("video has %d frames\n", nFrames);

    frames = cell(1, nFrames);
    framesDispPlot = cell(1, nFrames);
    for k = 1:nFrames
        fprintf("processing frame %d/%d...", k, nFrames);
        %get the frame, stabilize it, save it
        tgt = safeConvertFrame(readFrame(v));
        tgt = flat_field_correction(tgt, 10, 0, 'gaussianBlur');
        % tgt = adapthisteq(tgt, ...
        %     'ClipLimit', 0.05, ...
        %     'NumTiles', [16 16]); 
        [D, stabilized] = diffeomorphicDemon(tgt, refImg, tgt);
        frames{k} = stabilized;
        framesDispPlot{k} = D;
        disp("Done !");
    end

    arr2gif(frames, outGifPath);
    arr2gifPlots(framesDispPlot);
end

% helpers
function I = safeLoadGray(path)
    [img, map] = imread(path);

    if ~isempty(map)
        I = ind2gray(img, map);
    elseif ndims(img) == 3
        I = rgb2gray(img);
    else
        I = img;
    end

    if islogical(I), I = double(I); end
    if ~isfloat(I), I = im2double(I); end
end

function I = safeConvertFrame(frame)

    if size(frame, 3) == 3
        I = rgb2gray(frame);
    else
        I = frame;
    end

    I = im2double(I); % keep everything in [0,1]
end

function arr2gif(arr, filename)
    for k = 1:numel(arr)
        [A, map] = gray2ind(arr{k}, 256);
    
        if k == 1
            imwrite(A, map, filename, "gif", "LoopCount", Inf, "DelayTime", 0.05);
        else
            imwrite(A, map, filename, "gif", "WriteMode", "append", "DelayTime", 0.05);
        end
    end
end

function arr2gifPlots(arr)
    for k = 1:numel(arr)
        D = arr{k};
        Dc = complex(D(:,:,1), D(:,:,2));

        %map abs(Dc) to the min and max of the value of the image
        %maybe comapre to min and max of the video not the image to keep quantitative value idk
        absLog = log1p(abs(Dc));
        absLog = mat2gray(absLog);
        absLog = imadjust(absLog, [prctile(absLog(:), 1) prctile(absLog(:), 99);], [0 1]);
        
        % Convert to indexed
        [AAbs, mapAbs] = gray2ind(absLog, 256);
        [AAngle, mapAngle] = gray2ind(angle(Dc), 256);

        % Write to GIF
        if k == 1
            imwrite(AAbs, mapAbs, "abs.gif", "gif", "LoopCount", Inf, "DelayTime", 0.05);
            imwrite(AAngle, mapAngle, "angle.gif", "gif", "LoopCount", Inf, "DelayTime", 0.05);
        else
            imwrite(AAbs, mapAbs, "abs.gif", "gif", "WriteMode", "append", "DelayTime", 0.05);
            imwrite(AAngle, mapAngle, "angle.gif", "gif", "WriteMode", "append", "DelayTime", 0.05);
        end
    end

    analyzeFFT3D(arr);
end

function analyzeFFT3D(arr)
    % arr: cell array of displacement fields (MxNx2), one per frame

    nFrames = numel(arr);
    [H,W,~] = size(arr{1});
    stack = zeros(H, W, nFrames);
    for k = 1:nFrames
        D = arr{k};
        stack(:,:,k) = hypot(D(:,:,1), D(:,:,2));  % displacement magnitude
    end

    % F2Dc1 = fft(stack, [], 1);
    % F2Dc2 = fft(F2Dc1, [], 2);
    % F2Dc3 = fft(F2Dc2, [], 3);
    % 
    % s_omega = squeeze(mean(mean(abs(F2Dc3).^2, 1), 2));

    FD = fft2(BA(stack));
    FD = FD .* fftshift(disk);
    D_filt = ifft2(FD);

    figure;
    plot(log(s_omega));
    xlabel('Temporal frequency bin');
    ylabel('Power');
    title('Harmonic series of displacement magnitude');
end