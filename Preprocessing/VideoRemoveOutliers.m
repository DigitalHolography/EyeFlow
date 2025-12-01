function VideoRemoveOutliers(obj, params)
% Outlier Cleaning

% Compute the average profile over time (mean intensity per frame)
[numX, numY, ~] = size(obj.f_RMS);
diaphragmRadius = params.json.Mask.DiaphragmRadius;
maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
frame_means = squeeze(sum(obj.f_RMS .* maskDiaphragm, [1, 2]) / sum(maskDiaphragm, 'all')); % 1D array: numFrames x 1
mean_ = mean(frame_means);
std_ = std(frame_means);

% Detect outlier frames
outliers = (abs(frame_means - mean_) > 3 * std_);

kernel = [1 1 1 1 1 1 1 1]; % simple 1D smoothing
outliers_dilated = conv(double(outliers), kernel, 'same') > 0;

% Find the biggest contiguous part of outliers
[labeled, numRegions] = bwlabel(outliers_dilated);
if numRegions > 0
    regionSizes = histcounts(labeled, 1:numRegions+1);
    [~, largestIdx] = max(regionSizes);
    outlier_frames_mask = (labeled == largestIdx);
else
    outlier_frames_mask = false(size(outliers_dilated));
end

% If no outliers detected, return early
if ~any(outlier_frames_mask)
    return;
end

% Interpolate outlier frames for each video
obj.M0(:,:,outlier_frames_mask) = NaN;
obj.M1(:,:,outlier_frames_mask) = NaN;
obj.M2(:,:,outlier_frames_mask) = NaN;
obj.M0_ff(:,:,outlier_frames_mask) = NaN;
obj.f_RMS(:,:,outlier_frames_mask) = NaN;
obj.f_AVG(:,:,outlier_frames_mask) = NaN;

end
