function preMaskArtery = preMaskArtery(video, maskVesselness)
% video: 3D matrix (H x W x T)
% maskVesselness: 2D binary mask (H x W)
% preMaskArtery: 2D mask containing two branches with highest correlation

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% -------------------------------
% Step 1: Separate mask into branches
[label, n] = labelVesselBranches(maskVesselness, true(size(maskVesselness)), ToolBox.Cache.list.xy_barycenter);
imwrite(uint16(label), 'artery_16_PreMask_labels.png');

numBranches = max(label(:));
numFrames = size(video, 3);

signals = zeros(numBranches, numFrames);

% -------------------------------
% Step 2: Compute average signal of video for each branch
for i = 1:numBranches
    branchMask = (label == i); % logical mask of branch i
    % Extract all pixels for branch i across frames
    branchPixels = reshape(video(repmat(branchMask, [1 1 numFrames])), [], numFrames);
    % Average temporal signal
    signals(i, :) = mean(branchPixels, 1);
end

% -------------------------------
% Step 3: Normalize signals
signals_n = (signals - mean(signals, 2)) ./ std (signals, [], 2); % normalize each branch signal

s_idx = select_regular_peaks(signals_n,0.5,0.3);



% -------------------------------
% Step 4: Combine them into final mask
preMaskArtery = false(size(maskVesselness));

for i=1:length(s_idx)
    preMaskArtery = preMaskArtery | (label == i);
end

end
