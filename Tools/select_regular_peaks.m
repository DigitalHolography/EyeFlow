function s_idx = select_regular_peaks(signals_n, method, params)
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
    params.threshold (1, 1) double {mustBePositive}
    params.tolerance (1, 1) double {mustBePositive, mustBeLessThan(params.tolerance, 1)}
end

ToolBox = getGlobalToolBox;
dt = ToolBox.stride / ToolBox.fs / 1000;
fs = 1 / dt;

% Compute gradient of each signal
gradient_n = gradient(signals_n);

switch method
    case "minmax"
        s_idx = select_minmax(signals_n, gradient_n, fs, dt);
    case "regular"
        s_idx = select_regular(gradient_n, params.threshold, params.tolerance);
    case "kmeans_cosine"
        s_idx = select_kmeans(signals_n,'cosine');
end

end

function s_idx = select_minmax(signals_n, gradient_n, fs, dt)
% Select signals based on presence of regular min/max peaks in derivative

[numBranches, numFrames] = size(signals_n);

% Average normalized signal across all branches
avgSignal = mean(signals_n, 1);

% Compute FFT
Y = fft(avgSignal);
P2 = abs(Y / numFrames);
P1 = P2(1:floor(numFrames / 2) + 1);
P1(2:end - 1) = 2 * P1(2:end - 1);

% Frequency vector
f = fs * (0:(numFrames / 2)) / numFrames;

% Find dominant frequency in physiological range (e.g. 0.5 - 5 Hz)
f_range = (f > 0.5 & f < 5); % 30 - 300 bpm
[~, idx] = max(P1(f_range));
f0 = f(f_range);
f0 = f0(idx);
idx0 = round(f0 / dt);

% Select signals based on presence of peaks in derivative at intervals ~1/f0
s_idx = zeros(1, numBranches);

for i = 1:numBranches
    [peaks, locs] = findpeaks(abs(gradient_n(i, :)), 'MinPeakDistance', 0.6 * idx0);
    [~, loc_ind] = max(peaks);

    if gradient_n(i, locs(loc_ind)) > 0
        s_idx(i) = 1;
    else
        s_idx(i) = 0;
    end

end

end

function s_idx = select_regular(gradient_n, threshold, tolerance)

% Select signals based on regularity of peaks in derivative

numBranches = size(gradient_n, 1);
s_idx = zeros(1, numBranches);

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

    end

end

end


function s_idx = select_kmeans(signals_n, distance)
% Select signals based on kmeans clustering with given distance

[numBranches, numFrames] = size(signals_n);
ToolBox = getGlobalToolBox;
dt = ToolBox.stride / ToolBox.fs / 1000;
fs = 1 / dt;

[idxs, C] = kmedoids(signals_n,2,"Distance",distance) ;

s_idx = zeros(1, numBranches);

gradient_C = gradient(C,1);

% just to gets idx0 the period
% 
avgSignal = C(1,:);

% Compute FFT
Y = fft(avgSignal);
P2 = abs(Y / numFrames);
P1 = P2(1:floor(numFrames / 2) + 1);
P1(2:end - 1) = 2 * P1(2:end - 1);

% Frequency vector
f = fs * (0:(numFrames / 2)) / numFrames;

% Find dominant frequency in physiological range (e.g. 0.5 - 5 Hz)
f_range = (f > 0.5 & f < 5); % 30 - 300 bpm
[~, idx] = max(P1(f_range));
f0 = f(f_range);
f0 = f0(idx);
idx0 = round(f0 / dt);

%1 index
[peaks1, locs1] = findpeaks(abs(gradient_C(1, :)), 'MinPeakDistance', 0.6 * idx0);

[~, loc_ind1] = max(peaks1);

if gradient_C(1, locs1(loc_ind1)) > 0
    s_idx(idxs==1) = 1;
else
    s_idx(idxs==1) = 0;
end

%2 index
[peaks2, locs2] = findpeaks(abs(gradient_C(2, :)), 'MinPeakDistance', 0.6 * idx0);

[~, loc_ind2] = max(peaks2);

if gradient_C(2, locs2(loc_ind2)) > 0
    s_idx(idxs==2) = 1;
else
    s_idx(idxs==2) = 0;
end

end
