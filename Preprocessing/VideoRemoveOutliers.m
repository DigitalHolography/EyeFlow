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
outlier_frames_mask = isoutlier(frame_means);

% If no outliers detected, return early
if ~any(outlier_frames_mask)
    return;
end

% figure, plot(frame_means), hold on, scatter((1:length(frame_means)), frame_means .* outlier_frames_mask)

% Interpolate outlier frames for each video
obj.M0_data_video = interpolateOutlierFrames(obj.M0_data_video, outlier_frames_mask);
obj.M1_data_video = interpolateOutlierFrames(obj.M1_data_video, outlier_frames_mask);
obj.M2_data_video = interpolateOutlierFrames(obj.M2_data_video, outlier_frames_mask);
obj.M0_ff_video = interpolateOutlierFrames(obj.M0_ff_video, outlier_frames_mask);
obj.f_RMS_video = interpolateOutlierFrames(obj.f_RMS_video, outlier_frames_mask);
obj.f_AVG_video = interpolateOutlierFrames(obj.f_AVG_video, outlier_frames_mask);
end

function video_cleaned = interpolateOutlierFrames(video, outlier_frames_mask)
% Input:
%   video: 3D array (numX x numY x numFrames)
%   outlier_frames_mask: logical array (1 x numFrames), true for outlier frames

% Output:
%   video_cleaned: 3D array (numX x numY x numFrames), with outlier frames interpolated

video_cleaned = video; % Initialize with the original video

% Find the indices of outlier frames
outlier_indices = find(outlier_frames_mask)';

% For each outlier frame, interpolate linearly using neighboring frames
for idx = outlier_indices
    % Find the previous and next non-outlier frames
    prev_frame = find(~outlier_frames_mask(1:idx - 1), 1, 'last'); % Last non-outlier before idx
    next_frame = find(~outlier_frames_mask(idx + 1:end), 1, 'first') + idx; % First non-outlier after idx

    % Handle edge cases (e.g., first or last frame is an outlier)
    if isempty(prev_frame)
        prev_frame = next_frame; % If no previous frame, use the next frame
    end

    if isempty(next_frame)
        next_frame = prev_frame; % If no next frame, use the previous frame
    end

    if prev_frame == next_frame
        video_cleaned(:, :, idx) = video_cleaned(:, :, prev_frame);
    else
        % Linearly interpolate the outlier frame
        alpha = (idx - prev_frame) / (next_frame - prev_frame); % Interpolation weight
        video_cleaned(:, :, idx) = (1 - alpha) * video(:, :, prev_frame) + alpha * video(:, :, next_frame);
    end

end

end
