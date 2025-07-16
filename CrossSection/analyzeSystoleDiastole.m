function analyzeSystoleDiastole(sysIdx, diasIdx, v_RMS, locsLabel, maskLabel, ...
    numCircles, numBranches, ToolBox, initial, papillaDiameter, vesselName, numFrames)

% Analyze systole and diastole cross-sections
% Inputs:
%   sysIdx, diasIdx - indices for systole and diastole phases
%   v_RMS - velocity data
%   locsLabel, maskLabel - artery location and mask data
%   numCircles, numBranches - dimensions of artery data
%   ToolBox - toolbox parameters
%   initial - prefix for patch names
%   papillaDiameter - diameter parameter
%   vesselName - name of vessel for plots
%   numFrames - total number of frames
% Outputs:
%   diameter_sys, diameter_dias - diameter measurements for systole/diastole
%   Q_sys, Q_dias - flux measurements for systole/diastole
%   diameters_diff, diameters_se_diff - differences between systole/diastole

gap_threshold = 3;

% Find where the gaps between consecutive indices exceed the threshold
discontinuity_points_sys = find(diff(sysIdx) > gap_threshold);
discontinuity_points_dias = find(diff(diasIdx) > gap_threshold);

% Determine the start and end indices for each segment
starts_sys = [1, discontinuity_points_sys + 1];
ends_sys = [discontinuity_points_sys, length(sysIdx)];
starts_dias = [1, discontinuity_points_dias + 1];
ends_dias = [discontinuity_points_dias, length(diasIdx)];

% Create the cell array by splitting at discontinuities
systole_cell = arrayfun(@(s, e) sysIdx(s:e), starts_sys, ends_sys, 'UniformOutput', false);
diastole_cell = arrayfun(@(s, e) diasIdx(s:e), starts_dias, ends_dias, 'UniformOutput', false);

% Ensure we start with systole by removing any diastole segments that come before first systole
first_sys_idx = sysIdx(1); % First systole index
first_dias_idx = diasIdx(1); % First diastole index

% Remove diastole segments that start before first systole
if ~isempty(diastole_cell) && first_dias_idx < first_sys_idx
    % Find which diastole segments to keep (those starting after first_sys_idx)
    keep_dias = cellfun(@(x) x(1) >= first_sys_idx, diastole_cell);
    diastole_cell = diastole_cell(keep_dias);

    % Update starts and ends for diastole
    starts_dias = starts_dias(keep_dias);
    ends_dias = ends_dias(keep_dias);
end

% Initialize output variables
diameters_sys = cell(1, length(systole_cell));
diameters_se_sys = cell(1, length(systole_cell));
diameters_dias = cell(1, length(diastole_cell));
diameters_se_dias = cell(1, length(diastole_cell));

% Calculate mean times for each segment
t_systole = zeros(1, length(systole_cell));

for i = 1:length(systole_cell)
    t_systole(i) = round(mean(systole_cell{i}));
end

t_diastole = zeros(1, length(diastole_cell));

for i = 1:length(diastole_cell)
    t_diastole(i) = round(mean(diastole_cell{i}));
end

% Process systole segments
for i = 1:length(systole_cell)
    sysIdx_i = systole_cell{i};
    sys_v_RMS = v_RMS(:, :, sysIdx_i);

    % Initialize cells for arteries
    D_cell_sys = cell(numCircles, numBranches);
    D_se_cell_sys = cell(numCircles, numBranches);

    % Cross-Section Analysis of the arteries
    parfor c_idx = 1:numCircles

        for b_idx = 1:numBranches

            if ~isempty(locsLabel{c_idx, b_idx})
                patchName_sys = sprintf('%s%d_C%d_systole', initial, b_idx, c_idx);
                [results_sys] = crossSectionAnalysis2(ToolBox, locsLabel{c_idx, b_idx}, ...
                    maskLabel{c_idx, b_idx}, sys_v_RMS, patchName_sys, papillaDiameter);

                D_cell_sys{c_idx, b_idx} = results_sys.D;
                D_se_cell_sys{c_idx, b_idx} = results_sys.D_se;
            end

        end

    end

    diameters_sys{i} = D_cell_sys;
    diameters_se_sys{i} = D_se_cell_sys;
end

% Process diastole segments
for i = 1:length(diastole_cell)
    diasIdx_i = diastole_cell{i};
    dias_v_RMS = v_RMS(:, :, diasIdx_i);

    % Initialize cells for arteries
    D_cell_dias = cell(numCircles, numBranches);
    D_se_cell_dias = cell(numCircles, numBranches);

    % Cross-Section Analysis of the arteries
    parfor c_idx = 1:numCircles

        for b_idx = 1:numBranches

            if ~isempty(locsLabel{c_idx, b_idx})
                patchName_dias = sprintf('%s%d_C%d_diastole', initial, b_idx, c_idx);
                [results_dias] = crossSectionAnalysis2(ToolBox, locsLabel{c_idx, b_idx}, ...
                    maskLabel{c_idx, b_idx}, dias_v_RMS, patchName_dias, papillaDiameter);

                D_cell_dias{c_idx, b_idx} = results_dias.D;
                D_se_cell_dias{c_idx, b_idx} = results_dias.D_se;
            end

        end

    end

    diameters_dias{i} = D_cell_dias;
    diameters_se_dias{i} = D_se_cell_dias;
end

% Calculate differences between systole and diastole
diameter_sys_array = nan(length(systole_cell), numCircles, numBranches);
diameter_se_sys_array = nan(length(systole_cell), numCircles, numBranches);
diameter_dias_array = nan(length(diastole_cell), numCircles, numBranches);
diameter_se_dias_array = nan(length(diastole_cell), numCircles, numBranches);

for c_idx = 1:numCircles

    for b_idx = 1:numBranches

        for idx = 1:length(diastole_cell)

            try
                % Store individual measurements
                diameter_sys_array(idx, c_idx, b_idx) = diameters_sys{idx}{c_idx, b_idx};
                diameter_se_sys_array(idx, c_idx, b_idx) = diameters_se_sys{idx}{c_idx, b_idx};
                diameter_dias_array(idx, c_idx, b_idx) = diameters_dias{idx}{c_idx, b_idx};
                diameter_se_dias_array(idx, c_idx, b_idx) = diameters_se_dias{idx}{c_idx, b_idx};
            catch

            end

        end

    end

end

diameter_sys_mean = mean(sum(diameter_sys_array, 3, 'omitnan'), 2, 'omitnan');
diameter_dias_mean = mean(sum(diameter_dias_array, 3, 'omitnan'), 2, 'omitnan');
diameter_diff_mean = mean(diameter_sys_array - diameter_dias_array, [2 3], 'omitnan');

diameter_se_sys_mean = sqrt(sum(sqrt(sum(diameter_se_sys_array .^ 2, 3, 'omitnan')) .^ 2, 2, 'omitnan')) ./ numCircles;
diameter_se_dias_mean = sqrt(sum(sqrt(sum(diameter_se_dias_array .^ 2, 3, 'omitnan')) .^ 2, 2, 'omitnan')) ./ numCircles;
diameter_se_diff_mean = sqrt(sum(sqrt((sum(diameter_se_sys_array .^ 2, 3, 'omitnan') ./ numBranches) .^ 2 + ...
    (sum(diameter_se_dias_array .^ 2, 3, 'omitnan')) ./ numBranches) .^ 2, 2, 'omitnan')) ./ numCircles;

% Plot results
T = ToolBox.stride / ToolBox.fs / 1000;

figure, hold on
errorbar(t_systole * T, diameter_sys_mean, diameter_se_sys_mean, ...
    'k', 'LineStyle', 'none', ...
    'Marker', '^', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
errorbar(t_diastole * T, diameter_dias_mean, diameter_se_dias_mean, ...
    'k', 'LineStyle', 'none', ...
    'Marker', 'v', 'MarkerFaceColor', 'auto', 'LineWidth', 2)

yline(median(diameter_sys_mean), 'k--', sprintf('%.1f', median(diameter_sys_mean)), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2, 'FontSize', 14)
yline(median(diameter_dias_mean), 'k--', sprintf('%.1f', median(diameter_dias_mean)), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2, 'FontSize', 14)

box on,
set(gca, 'LineWidth', 2)
pbaspect([1.618 1 1])
xlabel("Time (s)")
ylabel("Lumen Diameter (µm)")
axis padded
xlim([0 numFrames * T])
exportgraphics(gca, fullfile(ToolBox.path_png, 'local', ...
    sprintf('%s_plot_diasys_diameter_%s.png', ToolBox.folder_name, vesselName)))

%%
figure, hold on
errorbar((t_systole + t_diastole) * T / 2, ...
    diameter_diff_mean, ...
    diameter_se_diff_mean, ...
    'k', 'LineStyle', 'none', 'Marker', 'diamond', 'MarkerFaceColor', 'auto', 'LineWidth', 2)

yline(median(diameter_diff_mean), 'k--', ...
    sprintf('%.1f', median(diameter_diff_mean)), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2, 'FontSize', 14)

box on,
set(gca, 'LineWidth', 2)
pbaspect([1.618 1 1])
xlabel("Time (s)")
ylabel("\Delta Lumen Diameter (µm)")
axis padded
axP = axis;
axis([0, numFrames * T, 0, axP(4)])
exportgraphics(gca, fullfile(ToolBox.path_png, 'local', ...
    sprintf('%s_plot_diasys_diameter_diff_%s.png', ToolBox.folder_name, vesselName)))

end
