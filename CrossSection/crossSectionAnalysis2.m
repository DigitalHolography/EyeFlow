function [results] = crossSectionAnalysis2(ToolBox, loc, mask, v_RMS, patchName, papillaDiameter)

% Perform cross-section analysis on blood vessels.
%
% Inputs:
%   ToolBox     - Struct, contains parameters and paths.
%   locs        - Nx2 array, locations of vessel centers.
%   mask        - 2D array, mask for the region of interest.
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
v_masked = squeeze(mean(v_RMS, 3)) .* mask;
v_masked(~mask) = NaN;

% Define sub-image dimensions
subImgHW = round(0.01 * size(v_masked, 1) * params.json.CrossSectionsAnalysis.ScaleFactorWidth);

% Define sub-image dimensions
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

% Crop and rotate sub-image
subImg = cropCircle(subImg);
[rotatedImg, tilt_angle] = rotateSubImage(subImg);
rotatedImg(rotatedImg <= 0) = NaN;
results.subImg_cell = rescale(rotatedImg);

% Compute the Vessel Cross Section
[D, D_se, A, A_se, c1, c2, rsquare] = computeVesselCrossSection(rotatedImg, patchName, ToolBox, papillaDiameter);
results.D = D;
results.D_se = D_se;
results.A = A;


[histo,edges] = histcounts(subImg(:), 6);
results.v_histo = {histo,edges};

% Generate figures
saveCrossSectionFigure(rotatedImg, c1, c2, ToolBox, patchName);

% Initialize rejected masks
rejected_masks = zeros(numX, numY, 3);

if rsquare < 0.6 || isnan(D)
    rejected_masks(:, :, 1) = mask; % Red
else
    rejected_masks(:, :, 2) = mask;% Green
end

% Compute blood volume rate and average velocity

for t = 1:numFrames
    xRange = max(round(-subImgHW / 2) + loc(1), 1):min(round(subImgHW / 2) + loc(1), numX);
    yRange = max(round(-subImgHW / 2) + loc(2), 1):min(round(subImgHW / 2) + loc(2), numY);
    tmp = v_RMS(:, :, t) .* mask;
    subFrame = tmp(yRange, xRange);

    if size(subFrame, 1) < length(xRange) || size(subFrame, 2) < length(yRange)
        xRange = round(-subImgHW / 2) + loc(1):round(subImgHW / 2) + loc(1);
        yRange = round(-subImgHW / 2) + loc(2):round(subImgHW / 2) + loc(2);
        tmp = NaN(length(xRange), length(yRange));
        tmp(1:size(subFrame, 1), 1:size(subFrame, 2)) = subFrame;
        subFrame = tmp;
        clear tmp
    end

    subFrame = cropCircle(subFrame);
    subFrame = imrotate(subFrame, tilt_angle, 'bilinear', 'crop');

    v_profile = mean(subFrame, 1, 'omitnan');
    v_cross = mean(subFrame(c1:c2, :), 2, 'omitnan');

    % Compute average velocity
    v = mean(v_profile(c1:c2));

    % Compute standard deviation of velocity
    v_se = std(v_cross);

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
    results.v_profiles{t} = mean(subFrame, 1);
    results.v_profiles_se{t} = std(subFrame, [], 1);
end

results.rejected_masks = rejected_masks;

close all;
end
