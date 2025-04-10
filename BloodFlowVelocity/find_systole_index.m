function [sys_index_list, fullPulse, sys_max_list, sys_min_list] = find_systole_index(video, maskArtery)
% This function extracts systole and diastole indices from a noisy video signal.
% Inputs:
%   - video: 3D video data (height x width x time)
%   - maskArtery: Binary mask indicating the artery region
% Outputs:
%   - sys_index_list: Indices of systole peaks
%   - fullPulse: The extracted pulse signal
%   - sys_max_list: Indices of local maxima within each cycle
%   - sys_min_list: Indices of local minima within each cycle

% Step 1: Extract the pulse signal from the video using the artery mask
fullPulse = squeeze(sum(video .* maskArtery, [1 2]) / nnz(maskArtery));

% Denoise
detrendedPulse = detrend(fullPulse);
noOutlier = filloutliers(detrendedPulse, "center");
smoothPulse = smoothdata(noOutlier);

% Step 2: Compute the derivative of the smoothed signal
[~, peaks_idx] = findpeaks(smoothPulse);
diff_signal = diff(smoothPulse);

% Step 3: Detect systole peaks using findpeaks
min_peak_height = max(diff_signal) * 0.5; % Adaptive threshold
min_peak_distance = floor(length(fullPulse) / length(peaks_idx)) * 0.6; % Minimum distance between peaks
[~, sys_index_list] = findpeaks(diff_signal, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

% Step 4: Validate and clean up the detected peaks
sys_index_list = validate_peaks(sys_index_list, 10);

mean_distance_between_peaks = mean(diff(sys_index_list));

% Step 5: Find local maxima and minima within each cycle
sys_max_list = zeros(1, numel(sys_index_list));
sys_min_list = zeros(1, numel(sys_index_list));

% Find the minimum before the first cycle
[~, amin] = min(fullPulse(1:sys_index_list(1)));

if abs(amin-sys_index_list(1)) < 0.05 * mean_distance_between_peaks
    sys_min_list(1) = NaN; % Ignore if too close to the first peak
else
    sys_min_list(1) = sys_index_list(1) + amin;
end


for i = 1:(numel(sys_index_list) - 1)
    % Find the maximum within the current cycle
    [~, amax] = max(fullPulse(sys_index_list(i):sys_index_list(i + 1) - 1));
    sys_max_list(i) = sys_index_list(i) + amax;
    
    % Find the minimum within the current cycle
    [~, amin] = min(fullPulse(sys_index_list(i):sys_index_list(i + 1) - 1));
    sys_min_list(i+1) = sys_index_list(i) + amin;
end

% Find the maximum after the end cycle
[~, amax] = max(fullPulse(sys_index_list(i+1):end));
sys_max_list(i+1) = sys_index_list(i+1) + amax;

sys_max_list = sys_max_list';
sys_min_list = sys_min_list';


% Step 6: Error handling
if isempty(sys_index_list)
    error('No systole peaks detected. Check signal quality or adjust parameters.');
end

end

%% **Validate Peaks (Removes peaks that are too close)**
function sys_index_list = validate_peaks(sys_index_list, min_distance)
i = 1;

while i < numel(sys_index_list)
    
    if sys_index_list(i + 1) - sys_index_list(i) < min_distance
        sys_index_list(i + 1) = [];
    else
        i = i + 1;
    end
    
end

end
