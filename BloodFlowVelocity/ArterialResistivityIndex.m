function [] = ArterialResistivityIndex(t, v_video, mask, sysIdx, diasIdx, name, folder)

ToolBox = getGlobalToolBox;
% Color Maps
v_video = v_video .* mask;

if contains(name, 'Vein')
    tmp = diasIdx;
    diasIdx = sysIdx;
    sysIdx = tmp; clear tmp
end

v_masked = v_video .* mask;
v_masked(~mask) = NaN;
signal = squeeze(sum(v_masked, [1, 2], 'omitnan') / nnz(mask))';
dsignal = squeeze(std(v_masked, [], [1, 2], 'omitnan'))';

% Compute mean and standard deviation of vSys and vDias
vSys = mean(v_video(:, :, sysIdx), 3);
vDias = mean(v_video(:, :, diasIdx), 3);

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

% Clip RI to [0, 1] and handle NaNs
RI = (vSys - vDias) ./ vSys;
RI(RI > 1) = 1;
RI(RI < 0) = 0;
RI(isnan(RI)) = 0;

% PI calculation
PI_mean = (vSys_mean - vDias_mean) / v_mean;

% Partial derivatives for PI uncertainty
dPI_dvSys = 1 / v_mean;
dPI_dvDias = -1 / v_mean;
dPI_dvMean = -(vSys_mean - vDias_mean) / (v_mean^2);
dPI = sqrt( (dPI_dvSys * vSys_std)^2 + (dPI_dvDias * vDias_std)^2 + (dPI_dvMean * std(signal))^2 );

% Clip PI to [0, -] and handle NaNs
PI = (vSys - vDias) ./ v_mean;
PI(PI < 0) = 0;
PI(isnan(PI)) = 0;

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

% Save txt
fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_main_outputs', '.txt')), 'a');

if strcmp(name, 'velocity')
    fprintf(fileID, 'Mean Velocity artery : %f (mm/s) \r\n', v_mean);
    fprintf(fileID, 'Max Velocity artery : %f (mm/s) \r\n', vSys_mean);
    fprintf(fileID, 'Min Velocity artery : %f (mm/s) \r\n', vDias_mean);
end

fprintf(fileID, 'Arterial Resistivity Index (%s) : %f  \r\n', name, RI_mean);
fprintf(fileID, 'Arterial Pulsatility Index (%s) : %f  \r\n', name, PI_mean);
fclose(fileID);

% Save json

if strcmp(name, 'velocity')
    ToolBox.outputs.MeanVelocityArtery = v_mean;
    ToolBox.outputs.MeanVelocityArtery_std = std(signal);
    ToolBox.outputs.SysVelocityArtery = vSys_mean;
    ToolBox.outputs.SysVelocityArtery_std = vSys_std;
    ToolBox.outputs.DiaVelocityArtery = vDias_mean;
    ToolBox.outputs.DiaVelocityArtery_std = vDias_std;
end

ToolBox.outputs.(sprintf('RI%s', name)) = RI_mean;
ToolBox.outputs.(sprintf('dRI%s', name)) = dRI;
ToolBox.outputs.(sprintf('PI%s', name)) = PI_mean;
ToolBox.outputs.(sprintf('dPI%s', name)) = dPI;

end
