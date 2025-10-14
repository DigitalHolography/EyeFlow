function [M0_Systole_img, M0_Diastole_img, M0_Systole_video, M0_Diastole_video, sysindexes, diasindexes] = compute_diasys(M0_video, maskArtery, export_folder)

arguments
    M0_video
    maskArtery
    export_folder = []
end

ToolBox = getGlobalToolBox;
[~, ~, numFrames] = size(M0_video);
fullTime = ToolBox.Cache.t;

cDark = [1 0 0];
cLight = [1 0.5 0.5];

[sys_index_list, fullPulse, ~, ~] = find_systole_index(M0_video, maskArtery);

fullPulse = fullPulse';

% Check if sys_index_list is empty
if isempty(sys_index_list)
    warning('sys_index_list is empty. Skipping systole and diastole computation.');

    [~, amin] = min(M0_video, [], 3);
    [~, amax] = max(M0_video, [], 3);

    M0_Systole_img = M0_video(amax, 3);
    M0_Diastole_img = M0_video(amin, 3);
    M0_Systole_video = M0_video;
    M0_Diastole_video = M0_video;
    return;
end

figure("Visible", "off", "Color", "w");
hold on
numSys = numel(sys_index_list); % number of systoles
fpCycle = round(numFrames / numSys); % Frames per cycle

X = [fullTime, flip(fullTime)];
Y = [fullPulse, zeros(1, length(fullPulse))];
fill(X, Y, cLight, 'EdgeColor', 'none')

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
        plot(fullTime(sys_range), fullPulse(sys_range), 'Color', cDark, 'LineWidth', 2)
        X = [fullTime(sys_range), flip(fullTime(sys_range))];
        Y = [fullPulse(sys_range), zeros(1, length(sys_range))];
        fill(X, Y, cDark, 'EdgeColor', 'none')
    catch
    end

end

plot(fullTime, fullPulse, 'k', 'LineWidth', 2)

% Ensure sysindexes and diaindexes are within the bounds of the video size
sysindexes = sysindexes(sysindexes >= 1 & sysindexes <= numFrames);

sysindexes = sort(unique(sysindexes));
diasindexes = setdiff(1:numFrames, sysindexes);

% Compute the mean images
M0_Systole_img = mean(M0_video(:, :, sysindexes), 3, 'omitnan');
M0_Diastole_img = mean(M0_video(:, :, diasindexes), 3, 'omitnan');

M0_Systole_video = M0_video(:, :, sysindexes);
M0_Diastole_video = M0_video(:, :, diasindexes);

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
