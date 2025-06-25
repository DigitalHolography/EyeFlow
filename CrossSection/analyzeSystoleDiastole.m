function [width_sys, width_dias, Q_sys, Q_dias, widths_diff, widths_se_diff] = ...
    analyzeSystoleDiastole(sysIdx, diasIdx, v_RMS, locsLabel, maskLabel, ...
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
%   width_sys, width_dias - width measurements for systole/diastole
%   Q_sys, Q_dias - flux measurements for systole/diastole
%   widths_diff, widths_se_diff - differences between systole/diastole

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
first_sys_idx = sysIdx(1);  % First systole index
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
width_sys = zeros(3, length(systole_cell));
width_dias = zeros(3, length(diastole_cell));
Q_sys = zeros(2, length(systole_cell));
Q_dias = zeros(2, length(systole_cell));

widths_sys = cell(1, length(systole_cell));
widths_se_sys = cell(1, length(systole_cell));
widths_dias = cell(1, length(diastole_cell));
widths_se_dias = cell(1, length(diastole_cell));

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
    A_cell_sys = cell(numCircles, numBranches);
    D_cell_sys = cell(numCircles, numBranches);
    D_se_cell_sys = cell(numCircles, numBranches);
    Q_cell_sys = cell(numCircles, numBranches);
    Q_se_cell_sys = cell(numCircles, numBranches);

    % Cross-Section Analysis of the arteries
    parfor c_idx = 1:numCircles

        for b_idx = 1:numBranches

            if ~isempty(locsLabel{c_idx, b_idx})
                patchName_sys = sprintf('%s%d_C%d_systole', initial, b_idx, c_idx);
                [results_sys] = crossSectionAnalysis2(ToolBox, locsLabel{c_idx, b_idx}, ...
                    maskLabel{c_idx, b_idx}, sys_v_RMS, patchName_sys, papillaDiameter);

                A_cell_sys{c_idx, b_idx} = results_sys.A;
                D_cell_sys{c_idx, b_idx} = results_sys.D;
                D_se_cell_sys{c_idx, b_idx} = results_sys.D_se;
                Q_cell_sys{c_idx, b_idx} = results_sys.Q;
                Q_se_cell_sys{c_idx, b_idx} = results_sys.Q_se;
            end

        end

    end

    [medianWidth, avgWidth, stdWidth] = widthHistogram(D_cell_sys, D_se_cell_sys, ...
        A_cell_sys, sprintf('%s_Systole_%d', vesselName, i));
    [radiusQ, radiusQSE] = averageRadiusFlux(Q_cell_sys, Q_se_cell_sys);

    widths_sys{i} = D_cell_sys;
    widths_se_sys{i} = D_se_cell_sys;
    Qs_sys{i} = Q_cell_sys;
    Qs_se_sys{i} = Q_se_cell_sys;

    width_sys(:, i) = [medianWidth; avgWidth; stdWidth];
    Q_sys(:, i) = [mean(radiusQ, 'all'); sqrt(sum(radiusQSE .* radiusQSE, 'all')) ./ numel(radiusQSE)];
end

% Process diastole segments
for i = 1:length(diastole_cell)
    diasIdx_i = diastole_cell{i};
    dias_v_RMS = v_RMS(:, :, diasIdx_i);

    % Initialize cells for arteries
    A_cell_dias = cell(numCircles, numBranches);
    D_cell_dias = cell(numCircles, numBranches);
    D_se_cell_dias = cell(numCircles, numBranches);
    Q_cell_dias = cell(numCircles, numBranches);
    Q_se_cell_dias = cell(numCircles, numBranches);

    % Cross-Section Analysis of the arteries
    parfor c_idx = 1:numCircles

        for b_idx = 1:numBranches

            if ~isempty(locsLabel{c_idx, b_idx})
                patchName_dias = sprintf('%s%d_C%d_diastole', initial, b_idx, c_idx);
                [results_dias] = crossSectionAnalysis2(ToolBox, locsLabel{c_idx, b_idx}, ...
                    maskLabel{c_idx, b_idx}, dias_v_RMS, patchName_dias, papillaDiameter);

                A_cell_dias{c_idx, b_idx} = results_dias.A;
                D_cell_dias{c_idx, b_idx} = results_dias.D;
                D_se_cell_dias{c_idx, b_idx} = results_dias.D_se;
                Q_cell_dias{c_idx, b_idx} = results_dias.Q;
                Q_se_cell_dias{c_idx, b_idx} = results_dias.Q_se;
            end

        end

    end

    [medianWidth, avgWidth, stdWidth] = widthHistogram(D_cell_dias, D_se_cell_dias, ...
        A_cell_dias, sprintf('%s_Diastole_%d', vesselName, i));
    [radiusQ, radiusQSE] = averageRadiusFlux(Q_cell_dias, Q_se_cell_dias);

    widths_dias{i} = D_cell_dias;
    widths_se_dias{i} = D_se_cell_dias;
    Qs_dias{i} = Q_cell_dias;
    Qs_se_dias{i} = Q_se_cell_dias;

    width_dias(:, i) = [medianWidth; avgWidth; stdWidth];
    Q_dias(:, i) = [mean(radiusQ, 'all'); sqrt(sum(radiusQSE .* radiusQSE, 'all')) ./ numel(radiusQSE)];
end

% Calculate differences between systole and diastole
widths_diff = nan(numCircles, numBranches);
widths_se_diff = nan(numCircles, numBranches);
Qs_diff = nan(length(diastole_cell), numCircles, numBranches);
Qs_se_diff = nan(numCircles, numBranches);

for c_idx = 1:numCircles

    for b_idx = 1:numBranches

        for idx = 1:length(diastole_cell)

            try
                widths_diff(c_idx, b_idx) = widths_sys{idx}{c_idx, b_idx} - widths_dias{idx}{c_idx, b_idx};
                widths_se_diff(c_idx, b_idx) = sqrt(widths_se_sys{idx}{c_idx, b_idx} .^ 2 + widths_se_dias{idx}{c_idx, b_idx} .^ 2);
                Qs_diff(idx, c_idx, b_idx) = mean(Qs_sys{idx}{c_idx, b_idx}) - mean(Qs_dias{idx}{c_idx, b_idx});
                % Q_se_diff(c_idx, b_idx) = sqrt(Qs_se_sys{idx}{c_idx, b_idx} .^ 2 + Qs_se_dias{idx}{c_idx, b_idx} .^ 2);
            end

        end

    end

end

width_diff = squeeze(mean(widths_diff, 1, 'omitnan'));
width_se_diff = squeeze(sqrt(sum(widths_se_diff .^ 2, 1, 'omitnan')) / numCircles);
Q_diff = squeeze(mean(sum(Qs_diff, 3), 2, 'omitnan'));

% Plot results
T = ToolBox.stride / ToolBox.fs / 1000;

figure, hold on
errorbar(t_systole * T, Q_sys(1, :), Q_sys(2, :), ...
    'k', 'LineStyle', 'none', ...
    'Marker', '^', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
errorbar(t_diastole * T, Q_dias(1, :), Q_dias(2, :), ...
    'k', 'LineStyle', 'none', ...
    'Marker', 'v', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
yline(mean(Q_sys(1, :)), 'k--', 'LineWidth', 2)
yline(mean(Q_dias(1, :)), 'k--', 'LineWidth', 2)

% errorbar((t_systole + t_diastole) * T / 2, ...
%     Q_sys(1, :) - Q_dias(1, :), ...
%     sqrt(Q_sys(2, :).^2 + Q_dias(2, :).^2), ...
%     'k', 'Marker', 'diamond', 'LineStyle', 'none', 'LineWidth', 2)
box on,
set(gca, 'LineWidth', 2)
pbaspect([1.618 1 1])
xlabel("Time (s)")
ylabel("Blood Volume Rate (µL/min)")
xlim([0 numFrames * T])
axis padded
exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', ...
    sprintf('%s_plot_diasys_flux_%s.png', ToolBox.folder_name, vesselName)))

figure, hold on
errorbar((t_systole + t_diastole) * T / 2, ...
    Q_diff, ...
    0, ...
    'k', 'Marker', 'diamond', 'LineStyle', 'none', 'LineWidth', 2)
box on,
set(gca, 'LineWidth', 2)
pbaspect([1.618 1 1])
xlabel("Time (s)")
ylabel("\Delta Blood Volume Rate (µL/min)")
xlim([0 numFrames * T])
axis padded
exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', ...
    sprintf('%s_plot_diasys_flux_diff_%s.png', ToolBox.folder_name, vesselName)))


figure, hold on
errorbar(t_systole * T, width_sys(1, :), width_sys(3, :), ...
    'k', 'LineStyle', 'none', ...
    'Marker', '^', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
errorbar(t_diastole * T, width_dias(1, :), width_dias(3, :), ...
    'k', 'LineStyle', 'none', ...
    'Marker', 'v', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
yline(median(width_sys(2, :)), 'k--', sprintf('%.1f', median(width_sys(2, :))), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2)
yline(median(width_dias(2, :)), 'k--', sprintf('%.1f', median(width_dias(2, :))), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2)
yline(median(width_sys(2, :) - width_dias(2, :)), 'k--', ...
    sprintf('%.1f', median(width_sys(2, :) - width_dias(2, :))), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2)
errorbar((t_systole + t_diastole) * T / 2, ...
    width_sys(2, :) - width_dias(2, :), ...
    sqrt(width_sys(3, :).^2 + width_dias(3, :).^2), ...
    'k', 'LineStyle', 'none', 'Marker', 'diamond', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
box on,
set(gca, 'LineWidth', 2)
pbaspect([1.618 1 1])
xlabel("Time (s)") 
ylabel("Lumen Diameter (µm)")
axis padded
xlim([0 numFrames * T])
exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', ...
    sprintf('%s_plot_diasys_width_%s.png', ToolBox.folder_name, vesselName)))

figure, hold on
errorbar((1:numBranches), width_diff, width_se_diff, 'k', ...
    'Marker', 'diamond', 'MarkerFaceColor', 'auto', 'LineWidth', 2)
box on,
set(gca, 'LineWidth', 2)
xlim([0 numBranches + 1])
pbaspect([1.618 1 1])
xlabel("Branch ID")
ylabel("\Delta Width (µm)")
exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', ...
    sprintf('%s_plot_diasys_width_diff_%s.png', ToolBox.folder_name, vesselName)))

end
