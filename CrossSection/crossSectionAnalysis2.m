function [results] = crossSectionAnalysis2(ToolBox, loc, ROI, xy_barycenter, v_RMS, patchName, papillaDiameter)

% Perform cross-section analysis on blood vessels.
%
% Inputs:
%   ToolBox     - Struct, contains parameters and paths.
%   locs        - Nx2 array, locations of vessel centers.
%   ROI         - 2D array, region of interest mask.
%   v_RMS       - 3D array, velocity data over time.
%   circleName  - String, name of the circle (for saving results).
%
% Outputs:
%   results     - Struct containing analysis results.

% Initialize parameters
params = ToolBox.getParams;
saveFigures = params.saveFigures;
[numX, numY, numFrames] = size(v_RMS);

% Initialize the results struct with preallocated fields.
results = struct();

% Compute mean velocity over time
v_masked = v_RMS;
v_masked(repmat(~ROI, [1, 1, size(v_RMS, 3)])) = NaN; % Apply mask to all slices

% Define sub-image dimensions
subImgHW = round(0.01 * size(v_masked, 1) * params.json.generateCrossSectionSignals.ScaleFactorWidth);

% Initialize results fields
xRange = max(round(-subImgHW / 2) + loc(1), 1):min(round(subImgHW / 2) + loc(1), numX);
yRange = max(round(-subImgHW / 2) + loc(2), 1):min(round(subImgHW / 2) + loc(2), numY);
subImg = v_masked(yRange, xRange, :);
subMask = ROI(yRange, xRange);

if size(subImg, 1) < length(xRange) || size(subImg, 2) < length(yRange)
    xRange = round(-subImgHW / 2) + loc(1):round(subImgHW / 2) + loc(1);
    yRange = round(-subImgHW / 2) + loc(2):round(subImgHW / 2) + loc(2);
    tmp = NaN(length(xRange), length(yRange), numFrames);
    tmp(1:size(subImg, 1), 1:size(subImg, 2), :) = subImg;
    subImg = tmp;
    clear tmp
end

subImgMean = squeeze(mean(subImg, 3, 'omitnan'));
subImgCropped = cropCircle(subImgMean);
[rotatedImg, tilt_angle] = rotateSubImage(subImgCropped);
subMask = imrotatecustom(subMask, tilt_angle);
rotatedImg(~subMask) = NaN;
results.subImg_cell = rescale(rotatedImg);

% UncroppedVersion
subImgUnCropped = squeeze(mean(v_RMS, 3));
subImgUnCropped = subImgUnCropped(yRange, xRange);
subImgUnCropped = imrotate(subImgUnCropped, tilt_angle, 'bilinear', 'crop');

% Compute the Vessel Cross Section
[D, D_SE, A, A_SE, c1, c2, rsquare] = computeVesselCrossSection(rotatedImg, params.px_size);
results.D = D;
results.D_SE = D_SE;
results.A = A;
results.A_SE = A_SE;

% Generate figures
if saveFigures
    saveCrossSectionFigure(subImgUnCropped, c1, c2, ToolBox, patchName);
end

% Initialize rejected masks
rejected_masks = zeros(numX, numY, 3);

if rsquare < 0.6 || isnan(D)
    rejected_masks(:, :, 1) = ROI; % Red
else
    rejected_masks(:, :, 2) = ROI; % Green
end

% Compute Flow Rate and average velocity

results.v_histo = cell(1, numFrames);

for t = 1:numFrames

    subFrame = subImg(:, :, t);
    subFrame = imrotatecustom(subFrame, tilt_angle);
    v_profile = mean(subFrame, 1, 'omitnan');
    v_cross = mean(subFrame(c1:c2, :), 2, 'omitnan');

    % Compute average velocity
    v = mean(v_profile(c1:c2), 'omitnan');

    [histo, ~] = histcounts(subFrame(~isnan(subFrame)), linspace(0, 60, 6)); % % HARD CODED
    results.v_histo{t} = histo;

    % Compute standard deviation of velocity
    v_SE = std(v_cross, 'omitnan');

    % Compute volumetric flow rate
    Q = v * A * 60; % microL/min

    % Uncertainty in volumetric flow rate
    if v ~= 0 && A ~= 0
        Q_SE = Q * sqrt((v_SE / v) ^ 2 + (A_SE / A) ^ 2 + (A_SE * v_SE / (A * v)) ^ 2);
    else
        Q_SE = 0; % Handle division by zero
    end

    % Handle NaN values
    if isnan(v)
        v = 0;
    end

    if isnan(Q)
        Q = 0;
    end

    if isnan(v_SE)
        v_SE = 0;
    end

    if isnan(Q_SE)
        Q_SE = 0;
    end

    % Store results
    results.v(t) = v;
    results.v_SE(t) = v_SE;
    results.Q(t) = Q;
    results.Q_SE(t) = Q_SE;
    results.v_profiles{t} = v_profile;
    results.v_profiles_SE{t} = std(subFrame, [], 1, 'omitnan');
end

results.rejected_masks = rejected_masks;

close all

end
