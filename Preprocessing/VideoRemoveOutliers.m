function VideoRemoveOutliers(obj, params)
% Outlier Cleaning

% Compute the average profile over time (mean intensity per frame)
[numX, numY, ~] = size(obj.f_RMS);
diaphragmRadius = params.json.Mask.DiaphragmRadius;
maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
frame_means = squeeze(sum(obj.f_RMS .* maskDiaphragm, [1, 2]) / sum(maskDiaphragm, 'all')); % 1D array: numFrames x 1

% Detect outlier frames
outlier_frames_mask = isoutlier(real(frame_means), "movmedian", 5, "ThresholdFactor", 2);

% figure, plot(frame_means), hold on, scatter((1:length(frame_means)), frame_means .* outlier_frames_mask)

% If no outliers detected, return early
if ~any(outlier_frames_mask)
    return;
end

% Interpolate outlier frames for each video
obj.M0 = interpolateOutlierFrames(obj.M0, outlier_frames_mask);
obj.M1 = interpolateOutlierFrames(obj.M1, outlier_frames_mask);
obj.M2 = interpolateOutlierFrames(obj.M2, outlier_frames_mask);
obj.M0_ff = interpolateOutlierFrames(obj.M0_ff, outlier_frames_mask);
obj.f_RMS = interpolateOutlierFrames(obj.f_RMS, outlier_frames_mask);
obj.f_AVG = interpolateOutlierFrames(obj.f_AVG, outlier_frames_mask);

end
