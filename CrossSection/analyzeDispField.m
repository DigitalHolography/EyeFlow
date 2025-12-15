function res = analyzeDispField(loc, ROI, xy_barycenter, dispField, patchName)
% Inputs:
%   ToolBox     - Struct, contains parameters and paths.
%   locs        - Nx2 array, locations of vessel centers.
%   ROI         - 2D array, region of interest mask.
%   dispField

    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;

    dispField = squeeze(hypot(dispField(:,:,1,:), dispField(:,:,2,:)));

    [numX, numY, numFrames] = size(dispField);

    % Compute mean velocity over time
    v_masked = dispField;
    v_masked(repmat(~ROI, [1, 1, size(patchName, 3)])) = NaN; % Apply mask to all slices

    subImgHW = round(0.01 * size(v_masked, 1) * params.json.generateCrossSectionSignals.ScaleFactorWidth);

    xRange = max(round(-subImgHW / 2) + loc(1), 1):min(round(subImgHW / 2) + loc(1), numX);
    yRange = max(round(-subImgHW / 2) + loc(2), 1):min(round(subImgHW / 2) + loc(2), numY);
    subImg = v_masked(yRange, xRange, :);
    subMask = ROI(yRange, xRange);

    subImgMean = squeeze(mean(subImg, 3, 'omitnan'));
    subImgCropped = cropCircle(subImgMean);
    [rotatedImg, tilt_angle] = rotateSubImage(subImgMean, subImgCropped, loc, xy_barycenter);
    % subMask = imrotatecustom(subMask, tilt_angle);
    % rotatedImg(~subMask) = NaN;

    for t = 1:numFrames
        subFrame = subImg(:, :, t);
        subFrame = imrotatecustom(subFrame, tilt_angle);
        res{t} = mean(subFrame, 1, 'omitnan');
    end


    % save history
    res_array = reshape(cell2mat(res), [17,256]);
    res_fft = fft(res_array, [], 2);
    % figure; plot(abs(res_fft(:, 2:end)))


end