function [preMaskArtery, preMaskVein] = preMaskArtery(video, maskVesselness)
% video: 3D matrix (H x W x T)
% maskVesselness: 2D binary mask (H x W)
% preMaskArtery: 2D mask containing two branches with highest correlation

ToolBox = getGlobalToolBox;
numFrames = size(video, 3);
fs = ToolBox.fs / ToolBox.stride * 1000; % Convert to seconds

% Step 1: Separate mask into branches
[label, ~] = labelVesselBranches(maskVesselness, true(size(maskVesselness)), ToolBox.Cache.xy_barycenter);
saveMaskImage(uint16(label), 'all_16_label_Vesselness.png', isStep = true);

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
[s_idx, locs_n] = select_regular_peaks(signals_n, 'minmax');

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

params = ToolBox.getParams;

% Find N for artery and N for vein

if params.saveFigures
    t = ToolBox.Cache.t;

    for N = 1:numBranches

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
        figure("Visible", 'off'); hold on;
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
        figPath = fullfile(ToolBox.path_png, 'mask', 'steps', sprintf('%s_%d_Peaks.png', name, N));
        exportgraphics(gcf, figPath, 'Resolution', 300);
    end

end

end
