function [preMaskArtery, preMaskVein] = preMaskArtery(video, maskVesselness)
% video: 3D matrix (H x W x T)
% maskVesselness: 2D binary mask (H x W)
% preMaskArtery: 2D mask containing two branches with highest correlation

ToolBox = getGlobalToolBox;
cmapArtery = ToolBox.Cache.cmapArtery;
numFrames = size(video, 3);
fs = ToolBox.fs / ToolBox.stride * 1000; % Convert to seconds

% Step 1: Separate mask into branches
[label, ~] = labelVesselBranches(maskVesselness, true(size(maskVesselness)), ToolBox.Cache.xy_barycenter);
saveMaskImage(uint16(label), 'artery_16_PreMask.png', isStep = true, cmap = cmapArtery);

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

% Step 3: Normalize signals
signals_n = (signals - mean(signals, 2)) ./ std (signals, [], 2); % normalize each branch signal
s_idx = select_regular_peaks(signals_n, 'minmax');

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
