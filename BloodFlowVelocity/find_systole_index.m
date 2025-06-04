function [sys_index_list, fullPulse, sys_max_list, sys_min_list] = find_systole_index(video, maskArtery, savepng)
% This function extracts systole and diastole indices from a noisy video signal.
% Inputs:
%   - video: 3D video data (height x width x time)
%   - maskArtery: Binary mask indicating the artery region
% Outputs:
%   - sys_index_list: Indices of systole peaks
%   - fullPulse: The extracted pulse signal
%   - sys_max_list: Indices of local maxima within each cycle
%   - sys_min_list: Indices of local minima within each cycle

% Input validation
if nargin < 3
    savepng = false;
end

% Step 1: Extract pulse signal
fullPulse = squeeze(sum(video .* maskArtery, [1 2], 'omitnan') / nnz(maskArtery));
fullPulseSmooth = smoothdata(filloutliers(fullPulse, "center", 'movmean', 5));

% Step 2: Compute derivative
diff_signal = gradient(fullPulseSmooth);

% Step 3: Detect peaks
min_peak_height = prctile(diff_signal, 80);
min_peak_distance = floor(length(fullPulse) / 10);
[~, sys_index_list] = findpeaks(diff_signal, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

% Step 4: Validate peaks
sys_index_list = validate_peaks(sys_index_list, 10);

% Step 5: Find local maxima and minima
num_peaks = numel(sys_index_list);
sys_max_list = zeros(num_peaks, 1);
sys_min_list = zeros(num_peaks, 1);

for i = 1:(numel(sys_index_list) - 1)
    L = sys_index_list(i + 1) - sys_index_list(i);
    D = round(L / 2);
    % Find the maximum within the current cycle
    [~, amax] = max(fullPulse(sys_index_list(i):sys_index_list(i) + D));
    sys_max_list(i) = sys_index_list(i) + amax - 1;

    % Find the minimum within the current cycle
    [~, amin] = min(fullPulse(sys_index_list(i) + D:sys_index_list(i + 1) - 1));
    sys_min_list(i + 1) = sys_index_list(i) + amin - 1 + D;
end

% Find the minimum before the first cycle
[~, amin] = min(fullPulse(1:sys_index_list(1)));
sys_min_list(1) = amin;

% Find the maximum after the end cycle
[~, amax] = max(fullPulse(sys_index_list(i + 1):end));
sys_max_list(i + 1) = sys_index_list(i + 1) + amax - 1;

sys_max_list = sys_max_list';
sys_min_list = sys_min_list';

% Step 6: Error handling
if isempty(sys_index_list)
    error('No systole peaks detected. Check signal quality or adjust parameters.');
end

% % FOR DEBUG
if savepng

    ToolBox = getGlobalToolBox();
    T = ToolBox.stride / ToolBox.fs / 1000;
    numFrames = size(video, 3);
    fullTime = linspace(0, numFrames * T, numFrames);

    figure(Visible = 'off');
    hold on
    plot(fullTime, diff_signal, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
    plot(fullTime, fullPulse, 'k-', 'LineWidth', 1.5);
    scatter((sys_max_list - 1) * T, fullPulse(sys_max_list), 'r', "filled")
    scatter((sys_min_list - 1) * T, fullPulse(sys_min_list), 'b', "filled")
    scatter((sys_index_list - 1) * T, fullPulse(sys_index_list), 'k', "filled")

    xline((sys_index_list - 1) * T, 'k--')
    xline((sys_min_list - 1) * T, 'b--')
    xline((sys_max_list - 1) * T, 'r--')
    hold off

    axis padded;
    axP = axis;
    axis tight;
    axT = axis;
    axis([axT(1), axT(2), axP(3), axP(4)]);
    box on
    set(gca, 'LineWidth', 2, 'PlotBoxAspect', [2.5 1 1])
    xlabel("Time (s)")
    ylabel("Velocity (mm/s)")

    if ~isfolder(fullfile(ToolBox.path_png, 'bloodFlowVelocity'))
        mkdir(fullfile(ToolBox.path_png, 'bloodFlowVelocity'))
    end

    exportgraphics(gca, fullfile(ToolBox.path_png, 'bloodFlowVelocity', sprintf("%s_%s", ToolBox.main_foldername, 'find_systoles_indices.png')))

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
