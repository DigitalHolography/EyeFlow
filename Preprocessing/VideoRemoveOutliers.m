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
N_ = params.json.Preprocess.RemoveOutliers.KeepBiggestContiguousStdDistance;
outliers = (abs(frame_means - mean_) > N_ * std_);

kernel = [1 1 1 1 1 1 1 1]; % simple 1D smoothing
outliers_dilated = conv(double(outliers), kernel, 'same') > 0;

% figure, plot(frame_means), hold on, yline(mean_ - N_ * std_), yline(mean_ + N_ * std_);

% If no outliers detected, return early
if ~any(outliers_dilated)
    return;
end

% If Keep biggest Contiguous Params is enabled, keep only the largest contiguous block of outliers
if params.json.Preprocess.RemoveOutliers.KeepBiggestContiguous
    % Find the biggest contiguous part of outliers
    
    i = find(outliers_dilated);
    if isempty(i)
        selected_range = true(size(outliers_dilated));
    else
        d = diff([0; i; numel(outliers_dilated) + 1]);
        [~, maxIdx] = max(d);
        if maxIdx == 1
            largest_contiguous_start = 1;
        else
            largest_contiguous_start = i(maxIdx - 1) + 1;
        end
        if maxIdx > numel(i)
            largest_contiguous_end = numel(outliers_dilated);
        else
            largest_contiguous_end = i(maxIdx) - 1;
        end
        selected_range = false(size(outliers_dilated));
        selected_range(largest_contiguous_start:largest_contiguous_end) = true;
    end

    obj.M0 = obj.M0(:, :, selected_range);
    obj.M1 = obj.M1(:, :, selected_range);
    obj.M2 = obj.M2(:, :, selected_range);
    obj.M0_ff = obj.M0_ff(:, :, selected_range);
    obj.f_RMS = obj.f_RMS(:, :, selected_range);
    obj.f_AVG = obj.f_AVG(:, :, selected_range);

elseif params.json.Preprocess.RemoveOutliers.InsertNans
    outlier_frames_mask = outliers_dilated;
    % Interpolate outlier frames for each video
    obj.M0(:, :, outlier_frames_mask) = NaN;
    obj.M1(:, :, outlier_frames_mask) = NaN;
    obj.M2(:, :, outlier_frames_mask) = NaN;
    obj.M0_ff(:, :, outlier_frames_mask) = NaN;
    obj.f_RMS(:, :, outlier_frames_mask) = NaN;
    obj.f_AVG(:, :, outlier_frames_mask) = NaN;
else
    outlier_frames_mask = outliers_dilated;

end

end
