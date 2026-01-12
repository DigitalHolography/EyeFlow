function stabilizeVideo(refImgPath, targetVideo, outGifPath)
refImg = safeLoadGray(refImgPath);
refImg = imresize(refImg, [512 512]);
refImg = flat_field_correction(refImg, 10, 0, 'gaussianBlur');

v = VideoReader(targetVideo);
nFrames = floor(v.Duration * v.FrameRate);
fprintf("video has %d frames\n", nFrames);

frames = cell(1, nFrames);
framesDispPlot = cell(1, nFrames);

for k = 1:nFrames
    fprintf("processing frame %d/%d...", k, nFrames);
    %get the frame, stabilize it, save it
    tgt = safeConvertFrame(readFrame(v));
    tgt = flat_field_correction_ef(tgt, 10, 0, 'gaussianBlur');
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

function FD = arr2gifPlots(arr)

for k = 1:numel(arr)
    D = arr{k};
    Dc = complex(D(:, :, 1), D(:, :, 2));

    %map abs(Dc) to the min and max of the value of the image
    %maybe comapre to min and max of the video not the image to keep quantitative value idk
    absLog = log1p(abs(Dc));
    absLog = mat2gray(absLog);
    absLog = imadjust(absLog, [prctile(absLog(:), 1) prctile(absLog(:), 99); ], [0 1]);

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

    % --- Low-pass filter ---
    [nx, ny] = size(absLog);
    FD = fft2(absLog);
    FD = fftshift(FD); % shift zero-freq to center

    [mask, ~] = diskMask(nx, ny, 0.02, 2); % keep high pass
    FD = FD .* mask; % apply filter

    % FD = ifftshift(FD);                        % shift back
    % D_filt = abs(ifft2(FD));                   % inverse FFT
    % framesLowPass{k} = D_filt;

end

%arr2gif(framesLowPass, "lowPassFilter.gif");

analyzeFFT3D(arr);
end

function analyzeFFT3D(arr)
% arr: cell array of displacement fields (MxNx2), one per frame

nFrames = numel(arr);
[H, W, ~] = size(arr{1});
stack = zeros(H, W, nFrames);

for k = 1:nFrames
    D = arr{k};
    stack(:, :, k) = complex(D(:, :, 1), D(:, :, 2));
    %absLog = log1p(abs(stack(:,:,k)));
    %absLog = mat2gray(absLog);
    %absLog = imadjust(absLog, [prctile(absLog(:), 1) prctile(absLog(:), 99);], [0 1]);
    %stack(:,:,k) = absLog .* exp(1j*angle(stack(:,:,k)));
end

FD3d = fftn(stack);

for k = 1:size(FD3d, 3)
    Dc = FD3d(:, :, k);

    %map abs(Dc) to the min and max of the value of the image
    %maybe comapre to min and max of the video not the image to keep quantitative value idk
    absLog = log1p(abs(Dc));
    absLog = mat2gray(absLog);
    absLog = imadjust(absLog, [prctile(absLog(:), 1) prctile(absLog(:), 99); ], [0 1]);

    % Convert to indexed
    [AAbs, mapAbs] = gray2ind(absLog, 256);
    [AAngle, mapAngle] = gray2ind(angle(Dc), 256);

    % Write to GIF
    if k == 1
        imwrite(AAbs, mapAbs, "absFD.gif", "gif", "LoopCount", Inf, "DelayTime", 0.05);
        imwrite(AAngle, mapAngle, "angleFD.gif", "gif", "LoopCount", Inf, "DelayTime", 0.05);
    else
        imwrite(AAbs, mapAbs, "absFD.gif", "gif", "WriteMode", "append", "DelayTime", 0.05);
        imwrite(AAngle, mapAngle, "angleFD.gif", "gif", "WriteMode", "append", "DelayTime", 0.05);
    end

    if (k == 5)
        figure("Visible", "off");
        imshow(AAbs);
    end

end

end
