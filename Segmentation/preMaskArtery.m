function [preMaskArtery, preMaskVein] = preMaskArtery(video, maskVesselness)
% video: 3D matrix (H x W x T)
% maskVesselness: 2D binary mask (H x W)
% preMaskArtery: 2D mask containing two branches with highest correlation

ToolBox = getGlobalToolBox;
numFrames = size(video, 3);
fs = ToolBox.fs / ToolBox.stride * 1000; % Convert to seconds
dt = 1 / fs;

% Step 1: Separate mask into branches
[label, n] = labelVesselBranches(maskVesselness, true(size(maskVesselness)), ...
    ToolBox.Cache.xy_barycenter, 'refine', false);
saveMaskImage(uint16(label), 'all_20_label_Vesselness.png', isStep = true, cmap = jet(n + 1));

% Step 2: Compute average signal of video for each branch

numBranches = max(label(:));
signals = zeros(numBranches, numFrames);

[b, a] = butter(4, 15 / (fs / 2), 'low');
movingMeanWindow = round(fs * 0.1); % ~0.1 s window; adjust as needed

parfor i = 1:numBranches
    branchMask = (label == i); % logical mask of branch i

    % Extract all pixels for branch i across frames
    branchPixels = reshape(video(repmat(branchMask, [1 1 numFrames])), [], numFrames);

    if isempty(branchPixels)
        continue
    end

    % Average all pixels in this branch
    sig = mean(branchPixels, 1);

    % Low-pass filter
    sig = filtfilt(b, a, sig);

    % Moving mean smoothing to further denoise the temporal signal
    if movingMeanWindow > 1
        sig = movmean(sig, movingMeanWindow);
    end

    signals(i, :) = sig;
end

% Step 3: Normalize signals
signals_n = (signals - mean(signals, 2)) ./ std (signals, [], 2); % normalize each branch signal
avgSignal = mean(signals_n, 1);

% Compute FFT
Y = fft(avgSignal);
P2 = abs(Y / numFrames);
P1 = P2(1:floor(numFrames / 2) + 1);
P1(2:end - 1) = 2 * P1(2:end - 1);

% Frequency vector
f = fft_freq_vector(fs, numFrames, true);

% Find dominant frequency in physiological range (e.g. 0.5 - 5 Hz)
f_range = (f > 0.5 & f < 5); % 30 - 300 bpm
[~, idx] = max(P1(f_range));
f0 = f(f_range);
f0 = f0(idx);
t0 = 1 / f0;
idx0 = round(t0 / dt);

[s_idx, locs_n] = select_regular_peaks(signals_n, 'minmax', idx0);

% Step 3.5: Check periodicity (autocorrelation-based)
isPeriodic = false(numBranches, 1);

parfor i = 1:numBranches
    sig = signals_n(i, :);
    isPeriodic(i) = check_validity(sig, fs, f0);
end

% Step 4: Combine them into final mask
preMaskArtery = false(size(maskVesselness));
preMaskVein = false(size(maskVesselness));

parfor i = 1:numBranches

    if ~isPeriodic(i)
        continue; % skip non-periodic branches
    end

    if s_idx(i) == 1
        preMaskArtery = preMaskArtery | (label == i);
    elseif s_idx(i) == 0
        preMaskVein = preMaskVein | (label == i);
    end

end

params = ToolBox.getParams;

if params.saveFigures
    t = ToolBox.Cache.t;
    path_png = ToolBox.path_png;
    parfor N = 1:numBranches

        if ~isPeriodic(N)
            continue; % skip non-periodic branches
        end

        if s_idx(N) == 1
            color = [1 0 0];
            name = 'artery';
        else
            color = [0 0 1];
            name = 'vein';
        end

        signals_N = signals_n(N, :);
        gradient_N = gradient(signals_N);
        locs_N = locs_n{N};

        f1 = figure("Visible", 'off'); hold on;
        plot(t, gradient_N, '-', ...
            'Color', color, 'LineWidth', 1.5);
        plot(t, signals_N, '--', ...
            'Color', color * 0.5 + 0.5, 'LineWidth', 1.5);
        scatter(t(locs_N), gradient_N(locs_N), ...
            'MarkerFaceColor', color, 'MarkerEdgeColor', color);

        axis padded
        axP = axis;
        axis tight
        axT = axis;
        axis([axT(1), axT(2), axP(3), axP(4)])
        box on
        set(gca, 'LineWidth', 2);
        set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
        figPath = fullfile(path_png, 'mask', 'steps', sprintf('%s_%d_Peaks.png', name, N));
        exportgraphics(gcf, figPath, 'Resolution', 300);
        close(f1)

    end

end

end

function isValid = check_validity(signal, fs, f0)
%CHECK_VALIDITY  Check if a temporal signal is periodic and not noise.
%
%   isValid = check_validity(signal, fs, f0)
%
%   INPUTS:
%       signal : 1 x N vector, normalized temporal signal of one branch
%       fs     : sampling frequency in Hz
%       f0     : reference (expected) fundamental frequency in Hz
%
%   OUTPUT:
%       isValid : true if the signal is periodic and matches f0,
%                 false otherwise.
%
%   The function uses the power spectrum of the signal to test:
%       - if it has sufficient total energy,
%       - if there is a strong dominant frequency,
%       - if that frequency lies close to f0.

% ---------------- Parameters ----------------
freqTolerance = 0.3; % Hz, tolerance around f0
dominanceThreshold = 3; % ratio of main peak power / average PSD
purityThreshold = 0.3;      % required purity (tune as needed)
freqRange = [0.5, 5]; % Hz, physiological range (30–300 bpm)

% ---------------- Preprocessing ----------------
signal = signal(:)' - mean(signal); % remove DC
numFrames = numel(signal);

% ---------------- Power Spectrum ----------------
Y = fft(signal);
P2 = abs(Y / numFrames) .^ 2; % power spectrum
P1 = P2(1:floor(numFrames / 2) + 1);
P1(2:end - 1) = 2 * P1(2:end - 1);
f = fs * (0:(numel(P1) - 1)) / numFrames;

% Restrict to physiological range
idxRange = (f >= freqRange(1) & f <= freqRange(2));
f_local = f(idxRange);
P_local = P1(idxRange);
P_local = P_local / sum(P_local); % normalize for purity calc

% ---- Dominant frequency ----
[~, idxMax] = max(P_local);
f_branch = f_local(idxMax);

% ---- Spectral Purity Metrics ----
% 1. Energy concentration near f0
band = (f_local > f_branch - 0.2 & f_local < f_branch + 0.2);
energyConcentration = sum(P_local(band));

% 2. Spectral entropy
spectralEntropy = -sum(P_local .* log(P_local + eps)) / log(numel(P_local));
purityEntropy = 1 - spectralEntropy;  % invert so 1 = pure, 0 = noisy

% 3. Combine into final purity score (weighted average)
purity = 0.7 * energyConcentration + 0.3 * purityEntropy;

% ---- Validity criteria ----
freqClose = abs(f_branch - f0) < freqTolerance;
dominance = P_local(idxMax) / mean(P_local);

isValid = freqClose && dominance > dominanceThreshold && purity > purityThreshold;

end
