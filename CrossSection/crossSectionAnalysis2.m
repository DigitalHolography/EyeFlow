function [results] = crossSectionAnalysis2(ToolBox, loc, ROI, v_RMS, patchName, papillaDiameter)

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
[numX, numY, numFrames] = size(v_RMS);

% Initialize the results struct with preallocated fields.
results = struct();

% Compute mean velocity over time
v_masked = squeeze(mean(v_RMS, 3)) .* ROI;
v_masked(~ROI) = NaN;

% Define sub-image dimensions
subImgHW = round(0.01 * size(v_masked, 1) * params.json.CrossSectionsAnalysis.ScaleFactorWidth);

% Initialize results fields
xRange = max(round(-subImgHW / 2) + loc(1), 1):min(round(subImgHW / 2) + loc(1), numX);
yRange = max(round(-subImgHW / 2) + loc(2), 1):min(round(subImgHW / 2) + loc(2), numY);
subImg = v_masked(yRange, xRange);

if size(subImg, 1) < length(xRange) || size(subImg, 2) < length(yRange)
    xRange = round(-subImgHW / 2) + loc(1):round(subImgHW / 2) + loc(1);
    yRange = round(-subImgHW / 2) + loc(2):round(subImgHW / 2) + loc(2);
    tmp = NaN(length(xRange), length(yRange));
    tmp(1:size(subImg, 1), 1:size(subImg, 2)) = subImg;
    subImg = tmp;
    clear tmp
end

% Interpolate the subImage two times
subImg = imresize(subImg, 2, 'bilinear');

% Crop and rotate sub-image
subImgCropped = cropCircle(subImg);
[rotatedImg, tilt_angle] = rotateSubImage(subImg,subImgCropped);
% rotatedImg(rotatedImg <= 0) = NaN;
results.subImg_cell = rescale(rotatedImg);

% Compute the Vessel Cross Section
[D, D_se, A, A_se, c1, c2, rsquare] = computeVesselCrossSection(rotatedImg, patchName, ToolBox, papillaDiameter);
results.D = D;
results.D_se = D_se;
results.A = A;

results.v_histo = cell(1, numFrames);

% Generate figures
subImgUnCropped = squeeze(mean(v_RMS, 3));
subImgUnCropped = subImgUnCropped(yRange, xRange);
subImgUnCropped = imresize(subImgUnCropped, 2, 'bilinear');
subImgUnCropped = imrotate(subImgUnCropped, tilt_angle, 'bilinear', 'crop');
saveCrossSectionFigure(subImgUnCropped, c1, c2, ToolBox, patchName);

% Initialize rejected masks
rejected_masks = zeros(numX, numY, 3);

if rsquare < 0.6 || isnan(D)
    rejected_masks(:, :, 1) = ROI; % Red
else
    rejected_masks(:, :, 2) = ROI; % Green
end

% Compute Volume Rate and average velocity

for t = 1:numFrames
    xRange = max(round(-subImgHW / 2) + loc(1), 1):min(round(subImgHW / 2) + loc(1), numX);
    yRange = max(round(-subImgHW / 2) + loc(2), 1):min(round(subImgHW / 2) + loc(2), numY);
    tmp = v_RMS(:, :, t) .* ROI;
    tmp(~ROI) = NaN;
    subFrame = tmp(yRange, xRange);

    if size(subFrame, 1) < length(xRange) || size(subFrame, 2) < length(yRange) % edge case (on the edges of the field)
        xRange = round(-subImgHW / 2) + loc(1):round(subImgHW / 2) + loc(1);
        yRange = round(-subImgHW / 2) + loc(2):round(subImgHW / 2) + loc(2);
        tmp = NaN(length(xRange), length(yRange));
        tmp(1:size(subFrame, 1), 1:size(subFrame, 2)) = subFrame;
        subFrame = tmp;
        clear tmp
    end

    subFrame = imresize(subFrame, 2, 'bilinear');

    subFrame = imrotatecustom(subFrame, tilt_angle);

    v_profile = mean(subFrame, 1, 'omitnan');
    v_cross = mean(subFrame(c1:c2, :), 2, 'omitnan');

    % Compute average velocity
    v = mean(subFrame(c1:c2, :),'all', 'omitnan');

    [histo, edges] = histcounts(subFrame(~isnan(subFrame)), linspace(0, 60, 6)); % % HARD CODED
    results.v_histo{t} = histo;

    % Compute standard deviation of velocity
    v_se = std(v_cross, 'omitnan');

    % Compute volumetric flow rate
    Q = v * A * 60; % microL/min

    % Uncertainty in volumetric flow rate
    if v ~= 0 && A ~= 0
        Q_se = Q * sqrt((v_se / v) ^ 2 + (A_se / A) ^ 2 + (A_se * v_se / (A * v)) ^ 2);
    else
        Q_se = 0; % Handle division by zero
    end

    % Handle NaN values
    if isnan(v)
        v = 0;
    end

    if isnan(Q)
        Q = 0;
    end

    if isnan(v_se)
        v_se = 0;
    end

    if isnan(Q_se)
        Q_se = 0;
    end

    % Store results
    results.v(t) = v;
    results.v_se(t) = v_se;
    results.Q(t) = Q;
    results.Q_se(t) = Q_se;
    results.v_profiles{t} = mean(subFrame, 1, 'omitnan');
    results.v_profiles_se{t} = std(subFrame, [], 1, 'omitnan');
end

results.rejected_masks = rejected_masks;

close all;
end
