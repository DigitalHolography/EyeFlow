function preMaskArtery = preMaskArtery(video, maskVesselness)
% video: 3D matrix (H x W x T)
% maskVesselness: 2D binary mask (H x W)
% preMaskArtery: 2D mask containing two branches with highest correlation

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% -------------------------------
% Step 1: Separate mask into branches
[~, label, ~] = getLongestArteryBranch(maskVesselness, ToolBox.Cache.list.xy_barycenter, 'all');
numBranches = max(label(:));
numFrames = size(video, 3);

% Preallocate
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
signalsNorm = zscore(signals, 0, 2); % normalize each branch signal

% Compute similarity matrix (dot products)
corrMatrix = signalsNorm * signalsNorm';

% Remove self-correlations
corrMatrix(1:numBranches + 1:end) = 0;

% Find the most correlated pair
[~, idx] = max(corrMatrix(:));
[branch1, branch2] = ind2sub(size(corrMatrix), idx);

% -------------------------------
% Step 4: Combine them into final mask
preMaskArtery = (label == branch1) | (label == branch2);
end
