function stabilizeVideo_test(refImgPath, targetVideo, outGifPath)
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
arr2gifPlots(framesDispPlot, v.framerate);
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

function arr2gifPlots(arr, framerate)

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

end

analyzeFFT3D(arr, framerate);
end

function analyzeFFT3D(arr, frameRate)
% arr: cell array of displacement fields (MxNx2), one per frame
% frameRate: The frame rate of the video in frames per second (Hz)

nFrames = numel(arr);

if nFrames < 2
    disp('Not enough frames for FFT analysis.');
    return;
end

[H, W, ~] = size(arr{1});
stack = zeros(H, W, nFrames);

% 1. Create the 3D stack of displacement magnitudes
for k = 1:nFrames
    D = arr{k};
    stack(:, :, k) = hypot(D(:, :, 1), D(:, :, 2)); % displacement magnitude
end

% 2. SUGGESTION: Remove the DC component (temporal mean) for each pixel
% This prevents a giant peak at 0 Hz in the final spectrum.
stack = stack - mean(stack, 3);

% 3. Perform the 3D FFT
% Applying fft along the 3rd dimension (time) is sufficient if we average later
F_time = fft(stack, [], 3);

% 4. Calculate the temporal power spectrum by averaging over space
% First, get power spectrum for each pixel: abs(F_time)^2
% Then, average across the spatial dimensions (1 and 2)
s_omega_unscaled = squeeze(mean(mean(abs(F_time) .^ 2, 1), 2));

% The spectrum is symmetric, so we only need the first half
powerSpectrum = s_omega_unscaled(1:floor(nFrames / 2) + 1);

% 5. Create the proper frequency axis in Hertz (Hz)
% This is the most critical step for interpretation.
f = frameRate * (0:(nFrames / 2)) / nFrames;

% 6. Find the peak frequency
[~, peakLocation] = max(powerSpectrum);
estimatedFreq = f(peakLocation);

fprintf('----------------------------------\n');
fprintf('FFT Analysis Complete.\n');
fprintf("frameRate: %.2f\n", frameRate);
fprintf('Estimated Dominant Frequency: %.2f Hz\n', estimatedFreq);
fprintf('Which is %.1f Beats Per Minute (BPM)\n', estimatedFreq * 60);
fprintf('----------------------------------\n');

% 7. Plot the results with the correct frequency axis
figure;
plot(f, powerSpectrum);
hold on;
plot(estimatedFreq, powerSpectrum(peakLocation), 'rv', 'MarkerFaceColor', 'r');
xlabel('Temporal Frequency (Hz)');
ylabel('Power');
title('Temporal Power Spectrum of Displacement Magnitude');
legend('', sprintf('Peak at %.2f Hz', estimatedFreq));
grid on;
xlim([0 frameRate / 2]); % Show up to the Nyquist frequency
end
