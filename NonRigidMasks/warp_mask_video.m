function warp_mask_video(refImgPath, refMaskPath, targetVideoPath, outGifPath)

    refImg  = safeLoadGray(refImgPath);
    refMask = safeLoadGray(refMaskPath);

    v = VideoReader(targetVideoPath);
    nFrames = floor(v.Duration * v.FrameRate);
    fprintf("video has %d frames\n", nFrames);

    first = safeConvertFrame(readFrame(v));
    refImg  = imresize(refImg,  size(first));
    refMask = imresize(refMask, size(first));
    v.CurrentTime = 0;

    [opt, metric] = imregconfig("monomodal");
    opt.MaximumIterations = 50;

    for k = 1:nFrames
        tgt = safeConvertFrame(readFrame(v));

        tform = imregtform(refImg, tgt, "rigid", opt, metric);
        R = imref2d(size(tgt));

        refW  = imwarp(refImg,  tform, "linear",  "OutputView", R);
        maskW = imwarp(refMask, tform, "nearest", "OutputView", R);  % preserve labels

        [D, ~] = imregdemons(refW, tgt, [50 25 10], ...
            "AccumulatedFieldSmoothing", 1.0, ...
            "PyramidLevels", 3, ...
            "DisplayWaitbar", false);

        maskW = imwarp(maskW, D, "nearest");

        if ~isequal(size(maskW), size(tgt))
            maskW = imresize(maskW, size(tgt));
        end

        frameU8 = im2uint8(maskW);
        [A, map] = gray2ind(frameU8, 256);

        if k == 1
            imwrite(A, map, outGifPath, "gif", "LoopCount", Inf, "DelayTime", 0.05);
        else
            imwrite(A, map, outGifPath, "gif", "WriteMode", "append", "DelayTime", 0.05);
        end

        fprintf("processed frame %d / %d\n", k, nFrames);
    end

    fprintf("saved warped mask GIF to: %s\n", outGifPath);
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
    if ~isfloat(I),  I = im2double(I); end
end

function I = safeConvertFrame(frame)
    if size(frame,3) == 3
        I = rgb2gray(frame);
    else
        I = frame;
    end
    I = im2double(I);  % keep everything in [0,1]
end
