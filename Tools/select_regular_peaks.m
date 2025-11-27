function [s_idx, locs_n] = select_regular_peaks(signals_n, method, idx0, params)
%SELECT_REGULAR_PEAKS Select signals with regular derivative peaks
%
% Inputs:
%   signals_n: normalized signals (numBranches x numFrames) double
%   method: 'regular' or 'minmax' (string)
%   params: struct with fields:
%       threshold: minimum peak height for 'regular' method (double)
%       tolerance: maximum coefficient of variation for 'regular' method (double)
%
% Outputs:
%   s_idx: indices of selected signals (1 x numBranches) double

arguments
    signals_n (:, :) double
    method (1, 1) string {mustBeMember(method, ["regular", "minmax", "kmeans_cosine"])}
    idx0
    params.threshold (1, 1) double {mustBePositive}
    params.tolerance (1, 1) double {mustBePositive, mustBeLessThan(params.tolerance, 1)}
end

% Compute gradient of each signal
gradient_n = gradient(signals_n);

switch method
    case "minmax"
        [s_idx, locs_n] = select_minmax(signals_n, gradient_n, idx0);
    case "regular"
        [s_idx, locs_n] = select_regular(gradient_n, params.threshold, params.tolerance);
    case "kmeans_cosine"
        s_idx = select_kmeans(signals_n, 'cosine');
end

end

function [s_idx, locs_n] = select_minmax(signals_n, gradient_n, idx0)
% Select signals based on presence of regular min/max peaks in derivative

[numBranches, ~] = size(signals_n);

% Select signals based on presence of peaks in derivative at intervals ~1/f0
s_idx = zeros(1, numBranches);
locs_n = cell(1, numBranches);

for i = 1:numBranches
    [peaks, locs] = findpeaks(abs(gradient_n(i, :)), 'MinPeakDistance', 0.8 * idx0);
    peaks_v = gradient_n(i, locs);
    c = sum(peaks_v > 0);

    if c > length(peaks) / 2
        s_idx(i) = 1;
    else
        s_idx(i) = 0;
    end

    locs_n{i} = locs;

end

end

function [s_idx, locs_n] = select_regular(gradient_n, threshold, tolerance)

% Select signals based on regularity of peaks in derivative

numBranches = size(gradient_n, 1);
s_idx = zeros(1, numBranches);
locs_n = cell(1, numBranches);

for i = 1:numBranches
    d_sig = gradient_n(i, :);

    if max(d_sig) > threshold

        % Find positive peaks in derivative
        [~, locs] = findpeaks(d_sig, 'MinPeakHeight', threshold);

        if numel(locs) > 2
            % Compute intervals between successive peaks
            intervals = diff(locs);

            % Check regularity: coefficient of variation (std/mean)
            cv = std(intervals) / mean(intervals);

            if cv < tolerance
                s_idx(i) = 1;
            end

        end

        locs_n{i} = locs;

    end

end

end

function s_idx = select_kmeans(signals_n, distance, idx0)
% Select signals based on kmeans clustering with given distance

[numBranches, ~] = size(signals_n);

[idxs, C] = kmedoids(signals_n, 2, "Distance", distance);

s_idx = zeros(1, numBranches);

gradient_C = gradient(C, 1);

%1 index
[peaks1, locs1] = findpeaks(abs(gradient_C(1, :)), 'MinPeakDistance', 0.6 * idx0);

[~, loc_ind1] = max(peaks1);

if gradient_C(1, locs1(loc_ind1)) > 0
    s_idx(idxs == 1) = 1;
else
    s_idx(idxs == 1) = 0;
end

%2 index
[peaks2, locs2] = findpeaks(abs(gradient_C(2, :)), 'MinPeakDistance', 0.6 * idx0);

[~, loc_ind2] = max(peaks2);

if gradient_C(2, locs2(loc_ind2)) > 0
    s_idx(idxs == 2) = 1;
else
    s_idx(idxs == 2) = 0;
end

end
