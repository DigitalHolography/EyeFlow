function [M0_Systole_img, M0_Diastole_img, M0_Systole_video, M0_Diastole_video, sysindexes, diasindexes] = compute_diasys(M0_ff_video, maskArtery, export_folder)

arguments
    M0_ff_video
    maskArtery
    export_folder = ''
end

ToolBox = getGlobalToolBox;
[~, ~, numFrames] = size(M0_ff_video);
fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

[sys_index_list, fullPulse, sys_max_list, sys_min_list] = find_systole_index(M0_ff_video, maskArtery);

fullPulse = fullPulse';

% Check if sys_index_list is empty
if isempty(sys_index_list)
    warning('sys_index_list is empty. Skipping systole and diastole computation.');

    [~, amin] = min(M0_ff_video, [], 3);
    [~, amax] = max(M0_ff_video, [], 3);

    M0_Systole_img = M0_ff_video(amax, 3);
    M0_Diastole_img = M0_ff_video(amin, 3);
    M0_Systole_video = M0_ff_video;
    M0_Diastole_video = M0_ff_video;
    return;
end

figAspect;
hold on

numSysMax = numel(sys_max_list); % number of systoles
numSysMin = numel(sys_min_list); % number of systoles
numSys = numel(sys_index_list); % number of systoles
fpCycleMax = round(numFrames / numSysMax); % Frames per cycle
fpCycleMin = round(numFrames / numSysMin); % Frames per cycle
fpCycle = round(numFrames / numSys); % Frames per cycle

sysindexes = [];
diasindexes = [];

for idx = 1:numSys
    xline(fullTime(sys_index_list(idx)), 'k--', 'LineWidth', 2)
    try
        % Calculate sysindexes and ensure the values stay within the valid range
        start_idx = sys_index_list(idx) + round(fpCycle * 0.05);
        end_idx = sys_index_list(idx) + round(fpCycle * 0.2);
        sys_range = start_idx:min(end_idx, numFrames);
        sysindexes = [sysindexes, sys_range];
        plot(fullTime(sys_range), fullPulse(sys_range), 'r-', 'LineWidth', 2)
        X = [fullTime(sys_range), flip(fullTime(sys_range))];
        Y = [fullPulse(sys_range), zeros(1, length(sys_range))];

        fill(X, Y, 'r', 'EdgeColor', 'none')
    catch
    end

    try
        start_idx = sys_index_list(idx) - round(fpCycle * 0.2);
        end_idx = sys_index_list(idx) - round(fpCycle * 0.05);
        dias_range = max(start_idx, 1):min(end_idx, numFrames);
        diasindexes = [diasindexes, dias_range];
        X = [fullTime(dias_range), flip(fullTime(dias_range))];
        Y = [fullPulse(dias_range), zeros(1, length(dias_range))];

        fill(X, Y, 'b', 'EdgeColor', 'none')
    catch
    end

end

for idx = 1:numSysMax

    try
        % Calculate sysindexes and ensure the values stay within the valid range
        start_idx = sys_max_list(idx);
        end_idx = sys_max_list(idx) + round(fpCycleMax * 0.15);
        sys_range = start_idx:min(end_idx, numFrames);
        sysindexes = [sysindexes, sys_range];

        X = [fullTime(sys_range), flip(fullTime(sys_range))];
        Y = [fullPulse(sys_range), zeros(1, length(sys_range))];

        fill(X, Y, 'r', 'EdgeColor', 'none')


    catch
    end

end

for idx = 1:numSysMin

    try
        % Calculate diaindexes and ensure the values stay within the valid range
        start_idx = sys_min_list(idx) - round(fpCycleMin * 0.15);
        end_idx = sys_min_list(idx);
        dias_range = max(start_idx, 1):min(end_idx, numFrames);
        diasindexes = [diasindexes, dias_range];
        plot(fullTime(dias_range), fullPulse(dias_range), 'b-', 'LineWidth', 2)

        X = [fullTime(dias_range), flip(fullTime(dias_range))];
        Y = [fullPulse(dias_range), zeros(1, length(dias_range))];

        fill(X, Y, 'b', 'EdgeColor', 'none')

    catch
    end

end

plot(fullTime, fullPulse, 'k', 'LineWidth', 2)

% Ensure sysindexes and diaindexes are within the bounds of the video size
sysindexes = sysindexes(sysindexes >= 1 & sysindexes <= numFrames);
diasindexes = diasindexes(diasindexes >= 1 & diasindexes <= numFrames);

sysindexes = sort(unique(sysindexes));
diasindexes = sort(unique(diasindexes));

% Compute the mean images
M0_Systole_img = mean(M0_ff_video(:, :, sysindexes), 3);
M0_Diastole_img = mean(M0_ff_video(:, :, diasindexes), 3);

M0_Systole_video = M0_ff_video(:, :, sysindexes);
M0_Diastole_video = M0_ff_video(:, :, diasindexes);

% Adjust axes
axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4) * 1.07])
title('diastole and systole')

xlabel('Time (s)')
pbaspect([1.618 1 1])

if strcmp(export_folder, 'mask')

    if isfolder(fullfile(ToolBox.path_png, 'mask'))
        ylabel('Power Doppler (a.u.)')
        exportgraphics(gca, fullfile(ToolBox.path_png, export_folder, 'steps', sprintf('%s_vessel_20_plot_diasys.png', ToolBox.main_foldername)))
    end

elseif strcmp(export_folder, 'bloodFlowVelocity')

    if isfolder(fullfile(ToolBox.path_png, 'bloodFlowVelocity'))
        ylabel('Velocity (mm/s)')
        exportgraphics(gca, fullfile(ToolBox.path_png, export_folder, sprintf('%s_diasysIdx.png', ToolBox.main_foldername)))
    end

end

end
