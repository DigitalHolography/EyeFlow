function [] = ArterialResistivityIndex(signal, systolesIndexes, name, folder, signal_se)

arguments
    signal
    systolesIndexes
    name
    folder
    signal_se = []
end

ToolBox = getGlobalToolBox;
numInterp = 60;
numFrames = length(signal);
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

if contains(name, 'velocity')
    unit = 'mm/s';
    y_label = 'Velocity (mm/s)';
else
    unit = 'µL/min';
    y_label = 'Blood Volume Rate (µL/min)';
end

% Color Maps

signal = double(signal);

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

% Formatting
if ~isempty(signal_se)
    title(sprintf('%s - RI: %.2f ± %.2f', name, RI, RI_se));
else
    title(sprintf('%s - RI: %.2f', name, RI));
end

xlabel('Time (s)');
ylabel(y_label);

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

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, folder, sprintf("%s_RI_%s.png", ToolBox.folder_name, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, folder, sprintf("%s_RI_%s.eps", ToolBox.folder_name, name)));
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

% Formatting
if ~isempty(signal_se)
    title(sprintf('%s - PI: %.2f ± %.2f', name, PI, PI_se));
else
    title(sprintf('%s - PI: %.2f', name, PI));
end

xlabel('Time (s)');
ylabel(y_label);

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), - 2, axP(4) + 2])

box on
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, folder, sprintf("%s_PI_%s.png", ToolBox.folder_name, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, folder, sprintf("%s_PI_%s.eps", ToolBox.folder_name, name)));
close;

if contains(name, 'Artery')
    VesselName = 'arterial';
else
    VesselName = 'venous';
end

% Save json

if contains(name, 'velocity')
    ToolBox.outputs.velocity.(sprintf('mean_%s', VesselName)) = round(v_mean, 2);
    ToolBox.outputs.velocity.(sprintf('systolic_%s', VesselName)) = round(vMax, 2);
    ToolBox.outputs.velocity.(sprintf('diastolic_%s', VesselName)) = round(vMin, 2);

    if ~isempty(signal_se)
        ToolBox.outputs.velocity.(sprintf('mean_%s_se', VesselName)) = round(std(interp_signal), 2);
        ToolBox.outputs.velocity.(sprintf('systolic_%s_se', VesselName)) = round(vMax_se, 2);
        ToolBox.outputs.velocity.(sprintf('diastolic_%s_se', VesselName)) = round(vMin_se, 2);
    end

    % New
    if contains(name, 'Vein')
        ToolBox.Outputs.add('VenousMeanVelocity', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('VenousMaximumVelocity', vMax, unit, vMax_se);
        ToolBox.Outputs.add('VenousMinimumVelocity', vMin, unit, vMin_se);
    elseif contains(name, 'Artery')
        ToolBox.Outputs.add('ArterialMeanVelocity', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('ArterialMinimumVelocity', vMin, unit, vMin_se);
        ToolBox.Outputs.add('ArterialMaximumVelocity', vMax, unit, vMax_se);
    end

else
    ToolBox.outputs.blood_volume_rate.(sprintf('mean_%s', VesselName)) = round(v_mean, 2);
    ToolBox.outputs.blood_volume_rate.(sprintf('systolic_%s', VesselName)) = round(vMax, 2);
    ToolBox.outputs.blood_volume_rate.(sprintf('diastolic_%s', VesselName)) = round(vMin, 2);

    if ~isempty(signal_se)
        ToolBox.outputs.blood_volume_rate.(sprintf('mean_%s_se', VesselName)) = round(std(interp_signal), 2);
        ToolBox.outputs.blood_volume_rate.(sprintf('systolic_%s_se', VesselName)) = round(vMax_se, 2);
        ToolBox.outputs.blood_volume_rate.(sprintf('diastolic_%s_se', VesselName)) = round(vMin_se, 2);
    end

    % New
    if contains(name, 'Vein')
        ToolBox.Outputs.add('VenousMeanVolumeRate', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('VenousMaximumVolumeRate', vMax, unit, vMax_se);
        ToolBox.Outputs.add('VenousMinimumVolumeRate', vMin, unit, vMin_se);
    elseif contains(name, 'Artery')
        ToolBox.Outputs.add('ArterialMeanVolumeRate', v_mean, unit, std(interp_signal));
        ToolBox.Outputs.add('ArterialMinimumVolumeRate', vMin, unit, vMin_se);
        ToolBox.Outputs.add('ArterialMaximumVolumeRate', vMax, unit, vMax_se);
    end

end

ToolBox.outputs.indices.(sprintf('%s_RI', name)) = round(RI, 2);
ToolBox.outputs.indices.(sprintf('%s_PI', name)) = round(PI, 2);
ToolBox.outputs.indices.(sprintf('%s_PR', name)) = round(PR, 2);

if ~isempty(signal_se)
    ToolBox.outputs.indices.(sprintf('%s_RI_se', name)) = round(RI_se, 2);
    ToolBox.outputs.indices.(sprintf('%s_PI_se', name)) = round(PI_se, 2);
    ToolBox.outputs.indices.(sprintf('%s_PR_se', name)) = round(PR_se, 2);
end

% New

if contains(name, 'velocity')

    if contains(name, 'Vein')
        ToolBox.Outputs.add('VenousResistivityIndexVelocity', RI, '', RI_se);
        ToolBox.Outputs.add('VenousPulsatilityIndexVelocity', PI, '', PI_se);
        ToolBox.Outputs.add('VenousMaxMinRatioVelocity', PR, '', PR_se);
    elseif contains(name, 'Artery')
        ToolBox.Outputs.add('ArterialResistivityIndexVelocity', RI, '', RI_se);
        ToolBox.Outputs.add('ArterialPulsatilityIndexVelocity', PI, '', PI_se);
        ToolBox.Outputs.add('ArterialMaxMinRatioVelocity', PR, '', PR_se);
    end

else

    if contains(name, 'Vein')
        ToolBox.Outputs.add('VenousResistivityIndexVolumeRate', RI, '', RI_se);
        ToolBox.Outputs.add('VenousPulsatilityIndexVolumeRate', PI, '', PI_se);
        ToolBox.Outputs.add('VenousMaxMinRatioVolumeRate', PR, '', PR_se);
    elseif contains(name, 'Artery')
        ToolBox.Outputs.add('ArterialResistivityIndexVolumeRate', RI, '', RI_se);
        ToolBox.Outputs.add('ArterialPulsatilityIndexVolumeRate', PI, '', PI_se);
        ToolBox.Outputs.add('ArterialMaxMinRatioVolumeRate', PR, '', PR_se);
    end

end

end
