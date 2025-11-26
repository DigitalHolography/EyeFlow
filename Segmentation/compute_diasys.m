function [M0_Systole_img, M0_Diastole_img, sysindexes, diasindexes] = compute_diasys(video, mask, export_folder)

arguments
    video
    mask
    export_folder = []
end

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;
[~, ~, numFrames] = size(video);
fullTime = ToolBox.Cache.t;

cDark = [1 0 0];
cLight = [1 0.5 0.5];

pulse_artery = squeeze(mean(video .* mask, [1 2], 'omitnan')) ./ nnz(sum(mask, [1 2]));

% Filter pulse_artery to remove high frequency noise
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
[b, a] = butter(4, 15 / (fs / 2), 'low');
pulse_artery = filtfilt(b, a, pulse_artery);

[sys_index_list, fullPulse, ~, ~] = find_systole_index(pulse_artery, 'savepng', false);

fullPulse = fullPulse';

% Check if sys_index_list is empty
if isempty(sys_index_list)
    warning('sys_index_list is empty. Skipping systole and diastole computation.');

    [~, amin] = min(video, [], 3);
    [~, amax] = max(video, [], 3);

    M0_Systole_img = video(amax, 3);
    M0_Diastole_img = video(amin, 3);
    return;
end

numSys = numel(sys_index_list); % number of systoles
fpCycle = round(numFrames / numSys); % Frames per cycle

if saveFigures
    figure("Visible", "off", "Color", "w");
    hold on

    X = [fullTime, flip(fullTime)];
    Y = [fullPulse, zeros(1, length(fullPulse))];
    fill(X, Y, cLight, 'EdgeColor', 'none')
end

% Find systole indexes and plot
sysindexes = [];

for idx = 1:numSys
    xline(fullTime(sys_index_list(idx)), 'k--', 'LineWidth', 2)

    try
        % Calculate sysindexes and ensure the values stay within the valid range
        start_idx = sys_index_list(idx);
        % end_idx = sys_index_list(idx) + round(fpCycle * 0.05);
        [~, end_idx] = max(fullPulse(start_idx:start_idx + round(fpCycle * 0.20)));
        end_idx = end_idx + start_idx;
        sys_range = start_idx:min(end_idx, numFrames);
        sysindexes = [sysindexes, sys_range];

        if saveFigures
            plot(fullTime(sys_range), fullPulse(sys_range), 'Color', cDark, 'LineWidth', 2)
            X = [fullTime(sys_range), flip(fullTime(sys_range))];
            Y = [fullPulse(sys_range), zeros(1, length(sys_range))];
            fill(X, Y, cDark, 'EdgeColor', 'none')
        end

    catch
    end

end

% Ensure sysindexes and diaindexes are within the bounds of the video size
sysindexes = sysindexes(sysindexes >= 1 & sysindexes <= numFrames);

sysindexes = sort(unique(sysindexes));
diasindexes = setdiff(1:numFrames, sysindexes);

% Compute the mean images
M0_Systole_img = mean(video(:, :, sysindexes), 3, 'omitnan');
M0_Diastole_img = mean(video(:, :, diasindexes), 3, 'omitnan');

if saveFigures
    plot(fullTime, fullPulse, 'k', 'LineWidth', 2)

    % Adjust axes
    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), 0, axP(4) * 1.07])

    xlabel('Time (s)')

    pbaspect([2.5 1 1]);
    box on
    set(gca, 'LineWidth', 2);
    ax = gca;
    ax.LineStyleOrderIndex = 1; % Reset if needed
    ax.SortMethod = 'depth'; % Try changing sorting method
    ax.Layer = 'top'; % This may help in some cases

    if isempty(export_folder)
        ylabel('Velocity (mm/s)')
        exportgraphics(gca, fullfile(ToolBox.path_png, sprintf('%s_diasysIdx.png', ToolBox.folder_name)))
    else
        ylabel('Power Doppler (a.u.)')
        exportgraphics(gca, fullfile(ToolBox.path_png, 'mask', 'steps', sprintf('%s_vessel_20_plot_diasys.png', ToolBox.folder_name)))
    end

end

end
