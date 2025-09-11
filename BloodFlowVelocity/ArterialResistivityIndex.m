function [] = ArterialResistivityIndex(signal, systolesIndexes, name, signal_se)

arguments
    signal
    systolesIndexes
    name
    signal_se = []
end

ToolBox = getGlobalToolBox;
numInterp = 60;
numFrames = length(signal);
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
fs = ToolBox.fs / ToolBox.stride * 1000; % Convert to seconds

if contains(name, 'v_')
    unit = 'mm/s';
    y_label = 'Velocity (mm/s)';
else
    unit = 'µL/min';
    y_label = 'Flow Rate (µL/min)';
end

% Color Maps

signal = double(signal);
[b, a] = butter(4, 15 / (fs / 2), 'low');
signal = filtfilt(b, a, signal);

[interp_signal, ~] = interpSignal(signal, systolesIndexes, numInterp);

% Compute mean and standard deviation of vSys and vDias
[~, amin] = min(interp_signal);
[~, amax] = max(interp_signal);

vMax = interp_signal(amax);
vMin = interp_signal(amin);

v_mean = mean(interp_signal);

% Compute RI and its uncertainty using error propagation
RI = (vMax - vMin) / vMax;

if ~isempty(signal_se)
    vMax_se = signal_se(amax);
    vMin_se = signal_se(amin);
    dRI_dvMax = vMin / (vMax ^ 2);
    dRI_dvMin = -1 / vMax;
    RI_se = sqrt((dRI_dvMax * vMax_se) ^ 2 + (dRI_dvMin * vMin_se) ^ 2);
else
    vMax_se = 0;
    vMin_se = 0;
    RI_se = 0;
end

% PI calculation
PI = (vMax - vMin) / v_mean;

if ~isempty(signal_se)
    dPI_dvMax = 1 / v_mean;
    dPI_dvMin = -1 / v_mean;
    dPI_dv_mean =- (vMax - vMin) / (v_mean ^ 2);
    PI_se = sqrt((dPI_dvMax * vMax_se) ^ 2 + (dPI_dvMin * vMin_se) ^ 2 + (dPI_dv_mean * std(interp_signal)) ^ 2);
else
    PI_se = 0;
end

% PR (Pulse Ratio) calculation
PR = vMax / vMin;

if ~isempty(signal_se)
    dPR_dvMax = 1 / vMin;
    dPR_dvMin = -vMax / (vMin ^ 2);
    PR_se = sqrt((dPR_dvMax * vMax_se) ^ 2 + (dPR_dvMin * vMin_se) ^ 2);
else
    PR_se = 0;
end

% RI Graph

figure('Visible', 'off');
hold on;

if ~isempty(signal_se)
    Color_std = [0.7, 0.7, 0.7];
    curve1 = signal + signal_se;
    curve2 = signal - signal_se;
    tmp_t = [t, fliplr(t)];
    inBetween = [curve1, fliplr(curve2)];
    fill(tmp_t, inBetween, Color_std);
    plot(t, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(t, curve2, "Color", Color_std, 'LineWidth', 2);
end

% Plot the velocity signal
plot(t, signal, '-k', 'LineWidth', 2);

% Add reference lines
yline(vMax, '--k', sprintf('%.1f %s', vMax, unit), 'LineWidth', 2, ...
    'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'top', 'FontWeight', 'bold');
yline(vMin, '--k', sprintf('%.1f %s', vMin, unit), 'LineWidth', 2, ...
    'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'bottom', 'FontWeight', 'bold');

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), - 2, axP(4) + 2])

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

% Add a text label with white background at the right edge of the plot
ax = gca;
xPos = ax.XLim(2) - ax.XLim(1); % Right edge of the plot
yLen = ax.YLim(2) - ax.YLim(1);

if ~isempty(signal_se)
    text(0.35 * xPos, 0.85 * yLen, sprintf("RI = %0.2f ± %0.2f", RI, RI_se), ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'Margin', 1); % Small padding
else
    text(0.4 * xPos, 0.85 * yLen, sprintf("RI = %0.2f", RI), ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'Margin', 1); % Small padding
end

xlabel('Time (s)');
ylabel(y_label);

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, sprintf("%s_RI_%s.png", ToolBox.folder_name, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, sprintf("%s_RI_%s.eps", ToolBox.folder_name, name)));
close;

% PI Graph

figure('Visible', 'off');
hold on;

if ~isempty(signal_se)
    Color_std = [0.7, 0.7, 0.7];
    curve1 = signal + signal_se;
    curve2 = signal - signal_se;
    tmp_t = [t, fliplr(t)];
    inBetween = [curve1, fliplr(curve2)];
    fill(tmp_t, inBetween, Color_std);
    plot(t, curve1, "Color", Color_std, 'LineWidth', 2);
    plot(t, curve2, "Color", Color_std, 'LineWidth', 2);
end

% Plot the velocity signal
plot(t, signal, '-k', 'LineWidth', 2);

% Add reference lines
yline(vMax, '--k', sprintf('%.1f %s', vMax, unit), 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'top', 'FontWeight', 'bold');
yline(v_mean, '--k', sprintf('%.1f %s', v_mean, unit), 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'FontWeight', 'bold');
yline(vMin, '--k', sprintf('%.1f %s', vMin, unit), 'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LabelVerticalAlignment', 'bottom', 'FontWeight', 'bold');

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), - 2, axP(4) + 2])

box on
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])

% Add a text label with white background at the right edge of the plot
ax = gca;
xPos = ax.XLim(2) - ax.XLim(1); % Right edge of the plot
yLen = ax.YLim(2) - ax.YLim(1);

if ~isempty(signal_se)
    text(0.35 * xPos, 0.85 * yLen, sprintf("PI = %.2f ± %.2f", PI, PI_se), ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'Margin', 1); % Small padding
else
    text(0.4 * xPos, 0.85 * yLen, sprintf("PI = %0.2f", PI), ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'Margin', 1); % Small padding
end

xlabel('Time (s)');
ylabel(y_label);


% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, sprintf("%s_PI_%s.png", ToolBox.folder_name, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, sprintf("%s_PI_%s.eps", ToolBox.folder_name, name)));
close;

if contains(name, 'artery')
    VesselName = 'arterial';
else
    VesselName = 'venous';
end

% Save json

if contains(name, 'v_')

    if contains(name, 'vein')
        ToolBox.Outputs.add('VenousMeanVelocity', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('VenousMaximumVelocity', vMax, unit, vMax_se);
        ToolBox.Outputs.add('VenousMinimumVelocity', vMin, unit, vMin_se);
    elseif contains(name, 'artery')
        ToolBox.Outputs.add('ArterialMeanVelocity', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('ArterialMinimumVelocity', vMin, unit, vMin_se);
        ToolBox.Outputs.add('ArterialMaximumVelocity', vMax, unit, vMax_se);
    end

else
    % New
    if contains(name, 'vein')
        ToolBox.Outputs.add('VenousMeanVolumeRate', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('VenousMaximumVolumeRate', vMax, unit, vMax_se);
        ToolBox.Outputs.add('VenousMinimumVolumeRate', vMin, unit, vMin_se);
    elseif contains(name, 'artery')
        ToolBox.Outputs.add('ArterialMeanVolumeRate', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('ArterialMinimumVolumeRate', vMin, unit, vMin_se);
        ToolBox.Outputs.add('ArterialMaximumVolumeRate', vMax, unit, vMax_se);
    end

end

if contains(name, 'v_')

    if contains(name, 'vein')
        ToolBox.Outputs.add('VenousResistivityIndexVelocity', RI, '', RI_se);
        ToolBox.Outputs.add('VenousPulsatilityIndexVelocity', PI, '', PI_se);
        ToolBox.Outputs.add('VenousMaxMinRatioVelocity', PR, '', PR_se);
    elseif contains(name, 'artery')
        ToolBox.Outputs.add('ArterialResistivityIndexVelocity', RI, '', RI_se);
        ToolBox.Outputs.add('ArterialPulsatilityIndexVelocity', PI, '', PI_se);
        ToolBox.Outputs.add('ArterialMaxMinRatioVelocity', PR, '', PR_se);
    end

else

    if contains(name, 'vein')
        ToolBox.Outputs.add('VenousResistivityIndexVolumeRate', RI, '', RI_se);
        ToolBox.Outputs.add('VenousPulsatilityIndexVolumeRate', PI, '', PI_se);
        ToolBox.Outputs.add('VenousMaxMinRatioVolumeRate', PR, '', PR_se);
    elseif contains(name, 'artery')
        ToolBox.Outputs.add('ArterialResistivityIndexVolumeRate', RI, '', RI_se);
        ToolBox.Outputs.add('ArterialPulsatilityIndexVolumeRate', PI, '', PI_se);
        ToolBox.Outputs.add('ArterialMaxMinRatioVolumeRate', PR, '', PR_se);
    end

end

end
