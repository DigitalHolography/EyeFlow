function obj = VideoRemoveOutliers(obj)
% Outlier Cleaning

params = Parameters_json(obj.directory, obj.param_name);

if ~params.json.Preprocess.RemoveOutliersFlag
    return
end

% Compute the average profile over time (mean intensity per frame)
[numX, numY, ~] = size(obj.f_RMS_video);
diaphragmRadius = params.json.Mask.DiaphragmRadius;
maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
frame_means = squeeze(sum(obj.f_RMS_video .* maskDiaphragm, [1, 2]) / sum(maskDiaphragm, 'all')); % 1D array: numFrames x 1

% Detect outlier frames
outlier_frames_mask = isoutlier(frame_means, "movmedian", 5, "ThresholdFactor", 2);

% figure, plot(frame_means), hold on, scatter((1:length(frame_means)), frame_means .* outlier_frames_mask)

% If no outliers detected, return early
if ~any(outlier_frames_mask)
    return;
end

% Interpolate outlier frames for each video
obj.M0_data_video = interpolateOutlierFrames(obj.M0_data_video, outlier_frames_mask);
obj.M1_data_video = interpolateOutlierFrames(obj.M1_data_video, outlier_frames_mask);
obj.M2_data_video = interpolateOutlierFrames(obj.M2_data_video, outlier_frames_mask);
obj.M0_ff_video = interpolateOutlierFrames(obj.M0_ff_video, outlier_frames_mask);
obj.f_RMS_video = interpolateOutlierFrames(obj.f_RMS_video, outlier_frames_mask);
obj.f_AVG_video = interpolateOutlierFrames(obj.f_AVG_video, outlier_frames_mask);

m = mean(obj.M0_ff_video, 'all');
s = std(obj.M0_ff_video, [], 'all');
obj.M0_ff_video(obj.M0_ff_video > m + 5 * s) = m + 5 * s;
obj.M0_ff_video(obj.M0_ff_video < m - 5 * s) = m - 5 * s;

end
