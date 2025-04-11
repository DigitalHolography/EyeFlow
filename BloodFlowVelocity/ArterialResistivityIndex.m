function [] = ArterialResistivityIndex(t, v_video, mask, sysIdx, diasIdx, name, folder)

ToolBox = getGlobalToolBox;
% Color Maps

if contains(name, 'Vein')
    tmp = diasIdx;
    diasIdx = sysIdx;
    sysIdx = tmp; clear tmp
end

if size(v_video, 3) ~= 1
    v_video = v_video .* mask;
    v_video(~mask) = NaN;
    signal = squeeze(sum(v_video, [1, 2], 'omitnan') / nnz(mask))';
    dsignal = squeeze(std(v_video, [], [1, 2], 'omitnan'))';

    % Compute mean and standard deviation of vSys and vDias
    vSys = mean(v_video(:, :, sysIdx), 3);
    vDias = mean(v_video(:, :, diasIdx), 3);
else
    signal = v_video;
    dsignal = mask;
end

vSys_frames = signal(sysIdx);
vDias_frames = signal(diasIdx);

vSys_mean = mean(vSys_frames);
vDias_mean = mean(vDias_frames);

vSys_std = std(vSys_frames);
vDias_std = std(vDias_frames);

v_mean = mean(signal);

% Compute RI and its uncertainty using error propagation
RI_mean = (vSys_mean - vDias_mean) / vSys_mean;

% Partial derivatives for error propagation
dRI_dvSys = vDias_mean / (vSys_mean^2);
dRI_dvDias = -1 / vSys_mean;

dRI = sqrt( (dRI_dvSys * vSys_std)^2 + (dRI_dvDias * vDias_std)^2 );

% PI calculation
PI_mean = (vSys_mean - vDias_mean) / v_mean;

% Partial derivatives for PI uncertainty
dPI_dvSys = 1 / v_mean;
dPI_dvDias = -1 / v_mean;
dPI_dvMean = -(vSys_mean - vDias_mean) / (v_mean^2);
dPI = sqrt( (dPI_dvSys * vSys_std)^2 + (dPI_dvDias * vDias_std)^2 + (dPI_dvMean * std(signal))^2 );

% PR (Pulse Ratio) calculation
PR_mean = vSys_mean / vDias_mean;

% Partial derivatives for PR uncertainty
dPR_dvSys = 1 / vDias_mean;
dPR_dvDias = -vSys_mean / (vDias_mean^2);
dPR = sqrt( (dPR_dvSys * vSys_std)^2 + (dPR_dvDias * vDias_std)^2 );

% RI Graph

figure('Visible', 'off');
hold on;

% Plot the uncertainty
Color_std = [0.7, 0.7, 0.7];
curve1 = signal + dsignal;
curve2 = signal - dsignal;
tmp_t = [t, fliplr(t)];
inBetween = [curve1, fliplr(curve2)];
fill(tmp_t, inBetween, Color_std);
plot(t, curve1, "Color", Color_std, 'LineWidth', 2);
plot(t, curve2, "Color", Color_std, 'LineWidth', 2);

% Plot the velocity signal
plot(t, signal, '-k', 'LineWidth', 2);

% Add reference lines
yline(vSys_mean, '--k', sprintf('%.1f mm/s', vSys_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
yline(vDias_mean, '--k', sprintf('%.1f mm/s', vDias_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

% Formatting
title(sprintf('%s - RI: %.2f ± %.2f', name, RI_mean, dRI));
xlabel('Time (s)');
ylabel('Velocity (mm/s)');

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4)])

fontsize(gca, 14, "points");
box on
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, folder, ...
    sprintf("%s_RI_%s.png", ToolBox.main_foldername, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, folder, ...
    sprintf("%s_RI_%s.eps", ToolBox.main_foldername, name)));
close;

% PI Graph

figure('Visible', 'off');
hold on;

% Plot the uncertainty
Color_std = [0.7, 0.7, 0.7];
curve1 = signal + dsignal;
curve2 = signal - dsignal;
tmp_t = [t, fliplr(t)];
inBetween = [curve1, fliplr(curve2)];
fill(tmp_t, inBetween, Color_std);
plot(t, curve1, "Color", Color_std, 'LineWidth', 2);
plot(t, curve2, "Color", Color_std, 'LineWidth', 2);

% Plot the velocity signal
plot(t, signal, '-k', 'LineWidth', 2);

% Add reference lines
yline(vSys_mean, '--k', sprintf('%.1f mm/s', vSys_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
yline(v_mean, '--k', sprintf('%.1f mm/s', v_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
yline(vDias_mean, '--k', sprintf('%.1f mm/s', vDias_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

% Formatting
title(sprintf('%s - PI: %.2f ± %.2f', name, PI_mean, dPI));
xlabel('Time (s)');
ylabel('Velocity (mm/s)');

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4)])

fontsize(gca, 14, "points");
box on
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, folder, ...
    sprintf("%s_PI_%s.png", ToolBox.main_foldername, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, folder, ...
    sprintf("%s_PI_%s.eps", ToolBox.main_foldername, name)));
close;

% Save image

if size(v_video, 3) > 1 % if given a video, output the image of RI / PI

    % Clip RI to [0, 1] and handle NaNs
    RI = (vSys - vDias) ./ vSys;
    RI(RI > 1) = 1;
    RI(RI < 0) = 0;
    RI(isnan(RI)) = 0;

    % Clip PI to [0, -] and handle NaNs
    PI = (vSys - vDias) ./ v_mean;
    PI(PI < 0) = 0;
    PI(isnan(PI)) = 0;

    % Generate colormap
    [cmap] = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);

    % Create RGB images for visualization

    % Display and save the RI image
    fig = figure("Visible", "off");
    imagesc(RI), axis image; axis off;
    colorbar, colormap(cmap)
    title(sprintf('RI %s = %0.2f', name, RI_mean));
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_RI_map_%s.png", ToolBox.main_foldername, name)));

    % Display and save the PI image
    f = figure("Visible", "off");
    imagesc(PI), axis image; axis off;
    colorbar, colormap(cmap)
    title(sprintf('PI %s = %0.2f', name, PI_mean));
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_PI_map_%s.png", ToolBox.main_foldername, name)));

    % Close figures
    close(f), close(fig);

else

end

if contains(name, 'Artery')
    VesselName = 'arterial';
else
    VesselName = 'venous';
end

% Save json

if contains(name, 'velocity')
    ToolBox.outputs.velocity.(sprintf('mean_%s', VesselName)) = round(v_mean, 2);
    ToolBox.outputs.velocity.(sprintf('mean_%s_se', VesselName)) = round(std(signal), 2);
    ToolBox.outputs.velocity.(sprintf('systolic_%s', VesselName))= round(vSys_mean, 2);
    ToolBox.outputs.velocity.(sprintf('systolic_%s_se', VesselName)) = round(vSys_std, 2);
    ToolBox.outputs.velocity.(sprintf('diastolic_%s', VesselName)) = round(vDias_mean, 2);
    ToolBox.outputs.velocity.(sprintf('diastolic_%s_se', VesselName)) = round(vDias_std, 2);

    % New
    if contains(name, 'Vein')
        ToolBox.Outputs.add('VenousMeanVelocity', v_mean, 'mm/s', std(signal));
        ToolBox.Outputs.add('VenousMaximumVelocity', vDias_mean, 'mm/s', vDias_std);
        ToolBox.Outputs.add('VenousMinimumVelocity', vSys_mean, 'mm/s', vSys_std);
    elseif contains(name, 'Artery')
        ToolBox.Outputs.add('ArterialMeanVelocity', v_mean, 'mm/s', std(signal));
        ToolBox.Outputs.add('ArterialMinimumVelocity', vDias_mean, 'mm/s', vDias_std);
        ToolBox.Outputs.add('ArterialMaximumVelocity', vSys_mean, 'mm/s', vSys_std);
    end
end

ToolBox.outputs.indices.(sprintf('%s_RI', name)) = round(RI_mean, 2);
ToolBox.outputs.indices.(sprintf('%s_RI_se', name)) = round(dRI, 2);
ToolBox.outputs.indices.(sprintf('%s_PI', name)) = round(PI_mean, 2);
ToolBox.outputs.indices.(sprintf('%s_PI_se', name)) = round(dPI, 2);
ToolBox.outputs.indices.(sprintf('%s_PR', name)) = round(PR_mean, 2);
ToolBox.outputs.indices.(sprintf('%s_PR_se', name)) = round(dPR, 2);

% New
if contains(name, 'Vein')
    ToolBox.Outputs.add('VenousResistivityIndexVelocity', RI_mean, '', dRI);
    ToolBox.Outputs.add('VenousPulsatilityIndexVelocity', PI_mean, '', dPI);
    ToolBox.Outputs.add('VenousMaxMinRatioVelocity', (vSys_mean / vDias_mean), '');
elseif contains(name, 'Artery')
    ToolBox.Outputs.add('ArterialResistivityIndexVelocity', RI_mean, '', dRI);
    ToolBox.Outputs.add('ArterialPulsatilityIndexVelocity', PI_mean, '', dPI);
    ToolBox.Outputs.add('ArterialMaxMinRatioVelocity', (vDias_mean / vSys_mean), '');
end
end
