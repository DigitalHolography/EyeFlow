function preMaskArtery = preMaskArtery(video, maskVesselness)
% video: 3D matrix (H x W x T)
% maskVesselness: 2D binary mask (H x W)
% preMaskArtery: 2D mask containing two branches with highest correlation

ToolBox = getGlobalToolBox;
% params = ToolBox.getParams;

numFrames = size(video, 3);
dt = ToolBox.stride / ToolBox.fs / 1000;
% t = linspace(0, numFrames * dt, numFrames);
fs = ToolBox.fs / ToolBox.stride * 1000; % Convert to seconds

% -------------------------------
% Step 1: Separate mask into branches
[label, n] = labelVesselBranches(maskVesselness, true(size(maskVesselness)), ToolBox.Cache.list.xy_barycenter);
imwrite(uint16(label), 'artery_16_PreMask_labels.png');

% -------------------------------
% Step 2: Compute average signal of video for each branch

numBranches = max(label(:));
signals = zeros(numBranches, numFrames);

for i = 1:numBranches
    branchMask = (label == i); % logical mask of branch i
    % Extract all pixels for branch i across frames
    branchPixels = reshape(video(repmat(branchMask, [1 1 numFrames])), [], numFrames);
    % Filter
    [b, a] = butter(4, 15 / (fs / 2), 'low');
    signals(i, :) = filtfilt(b, a, mean(branchPixels, 1));
end

% -------------------------------
% Step 3: Normalize signals
signals_n = (signals - mean(signals, 2)) ./ std (signals, [], 2); % normalize each branch signal
gradient_n = gradient(signals_n);

% Average normalized signal across all branches
avgSignal = mean(signals_n, 1);

% Compute FFT
N = length(avgSignal);
Y = fft(avgSignal);
P2 = abs(Y / N);
P1 = P2(1:floor(N / 2) + 1);
P1(2:end - 1) = 2 * P1(2:end - 1);

% Frequency vector
f = fs * (0:(N / 2)) / N;

% Find dominant frequency in physiological range (e.g. 0.5 - 5 Hz)
f_range = (f > 0.5 & f < 5); % 30 - 300 bpm
[~, idx] = max(P1(f_range));
f0 = f(f_range);
f0 = f0(idx);
idx0 = round(f0 / dt);

s_idx = zeros(1, n);

for i = 1:numBranches
    [peaks, locs] = findpeaks(abs(gradient_n(i, :)), 'MinPeakDistance', 0.6 * idx0);
    [~, loc_ind] = max(peaks);

    if gradient_n(i, locs(loc_ind)) > 0
        s_idx(i) = 1;
    else
        s_idx(i) = 0;
    end

end

s_idx = select_regular_peaks(signals_n, 0.5, 0.3);

% -------------------------------
% Step 4: Combine them into final mask
preMaskArtery = false(size(maskVesselness));

for i = 1:length(s_idx)

    if s_idx(i) == 1
        preMaskArtery = preMaskArtery | (label == i);
    end

end

preMaskVein = false(size(maskVesselness));

for i = 1:length(s_idx)

    if s_idx(i) == 0
        preMaskVein = preMaskVein | (label == i);
    end

end

end
