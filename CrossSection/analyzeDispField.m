function res = analyzeDispField(loc, ROI, xy_barycenter, dispField, ~)
% Inputs:
%   ToolBox     - Struct, contains parameters and paths.
%   locs        - Nx2 array, locations of vessel centers.
%   ROI         - 2D array, region of interest mask.
%   dispField

res = analyze_logic(loc, ROI, xy_barycenter, dispField);

% figure('Name', 'Line Profiles', 'Color', 'w');
%
% % Plot D
% subplot(3, 1, 1);
% plot(res.D); % Transpose to 256x33
% title('DispField');
% grid on;
%
% % Plot FFT
% subplot(3, 1, 2);
% plot(abs(res.fft_D));
% title('FFT');
% grid on;
%
% % Plot Derivative
% subplot(3, 1, 3);
% plot(res.dD_dt);
% title('dD/dt');
% grid on;

end

function res = analyze_logic(loc, ROI, xy_barycenter, dispField)
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% v_masked_resized = upscaleDispField(dispField);

[numX, numY, numFrames] = size(dispField);

% Compute mean velocity over time
v_masked = dispField;

% Apply the dilation
se = strel("disk", 5);
ROI_dilated = imdilate(ROI, se);
% v_masked(repmat(~ROI_dilated, [1, 1, numFrames])) = NaN; % Apply mask to all slices

subImgHW = round(0.01 * size(v_masked, 1) * params.json.generateCrossSectionSignals.ScaleFactorWidth);

xRange = max(round(-subImgHW / 2) + loc(1), 1):min(round(subImgHW / 2) + loc(1), numX);
yRange = max(round(-subImgHW / 2) + loc(2), 1):min(round(subImgHW / 2) + loc(2), numY);
subImg = dispField(yRange, xRange, :);
subMask = ROI_dilated(yRange, xRange);

% Apply the mask to only the relevant part
subImg(repmat(~subMask, [1, 1, numFrames])) = NaN;

if size(subImg, 1) < length(xRange) || size(subImg, 2) < length(yRange)
    xRange = round(-subImgHW / 2) + loc(1):round(subImgHW / 2) + loc(1);
    yRange = round(-subImgHW / 2) + loc(2):round(subImgHW / 2) + loc(2);
    tmp = NaN(length(xRange), length(yRange), numFrames);
    tmp(1:size(subImg, 1), 1:size(subImg, 2), :) = subImg;
    subImg = tmp;
    clear tmp
end

% subMask = ROI(yRange, xRange);

subImgMean = squeeze(mean(subImg, 3, 'omitnan'));
subImgCropped = cropCircle(subImgMean);

try
    [~, tilt_angle] = rotateSubImage(subImgMean, subImgCropped, loc, xy_barycenter);
catch ME
    printf("FAILED");
end

% subMask = imrotatecustom(subMask, tilt_angle);
% rotatedImg(~subMask) = NaN;

parfor t = 1:numFrames
    subFrame = subImg(:, :, t);
    subFrame = imrotatecustom(subFrame, tilt_angle);
    res_cell{t} = mean(subFrame, 1, 'omitnan');
end

% save history
res.D = reshape(cell2mat(res_cell), [size(xRange, 2), numFrames]);
res.fft_D = fft(res.D, [], 2);
% figure; plot(abs(res_fft(:, 2:end)))

[res.dD_dt, res.dD_dt_tmp] = gradient(res.D, 0.01, 1);
end

function res = handlePeaks(data)
res = [NaN, NaN, NaN, NaN];
[pks, locs] = findpeaks(data, "SortStr", "descend", "NPeaks", 2);

num_found = length(pks);

if num_found >= 1
    res(1) = locs(1);
    res(2) = pks(1);
end

if num_found >= 2
    res(3) = locs(2);
    res(4) = pks(2);
end

end
