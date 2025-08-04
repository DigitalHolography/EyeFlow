function [] = VenousResistivityIndex(t, v_video, mask, sysIdxList, name)

ToolBox = getGlobalToolBox;
numInterp = 60;
% Color Maps

if size(v_video, 3) ~= 1
    v_video = v_video .* mask;
    v_video(~mask) = NaN;
    signal = squeeze(sum(v_video, [1, 2], 'omitnan') / nnz(mask))';
    dsignal = squeeze(std(v_video, [], [1, 2], 'omitnan'))';
else
    signal = v_video;
    dsignal = mask;
end

try
    signal = double(wdenoise(signal, 4));
catch
    signal = double(signal);
end

[interp_signal, avgLength, interp_signal_ste] = interpSignal(signal, systolesIndexes, numInterp, dsignal);
interp_t = linspace(0, avgLength, numInterp);

vMax = max(interp_signal);
vMin = min(interp_signal);
vMax_frames_indxs = abs(interp_signal - vMax) / abs(vMax) < 0.10; % use value close to the min and max for calculating pulsatility
vMin_frames_indxs = abs(interp_signal - vMin) / abs(vMin) < 0.10;

vMax_mean = mean(interp_signal(vMax_frames_indxs));
vMin_mean = mean(interp_signal(vMin_frames_indxs));

vMax_std = std(interp_signal(vMax_frames_indxs));
vMin_std = std(interp_signal(vMin_frames_indxs));

v_mean = mean(signal);

% Compute RI and its uncertainty using error propagation
RI_mean = (vMax_mean - vMin_mean) / vMax_mean;

% Partial derivatives for error propagation
dRI_dvMax = vMin_mean / (vMax_mean ^ 2);
dRI_dvMin = -1 / vMax_mean;

dRI = sqrt((dRI_dvMax * vMax_std) ^ 2 + (dRI_dvMin * vMin_std) ^ 2);

% PI calculation
PI_mean = (vMax_mean - vMin_mean) / v_mean;

% Partial derivatives for PI uncertainty
dPI_dvMax = 1 / v_mean;
dPI_dvMin = -1 / v_mean;
dPI_dvMean =- (vMax_mean - vMin_mean) / (v_mean ^ 2);
dPI = sqrt((dPI_dvMax * vMax_std) ^ 2 + (dPI_dvMin * vMin_std) ^ 2 + (dPI_dvMean * std(signal)) ^ 2);

% PR (Pulse Ratio) calculation
PR_mean = vMax_mean / vMin_mean;

% Partial derivatives for PR uncertainty
dPR_dvMax = 1 / vMin_mean;
dPR_dvMin = -vMax_mean / (vMin_mean ^ 2);
dPR = sqrt((dPR_dvMax * vMax_std) ^ 2 + (dPR_dvMin * vMin_std) ^ 2);

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
yline(vMax_mean, '--k', sprintf('%.1f mm/s', vMax_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
yline(vMin_mean, '--k', sprintf('%.1f mm/s', vMin_mean), ...
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

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, ...
    sprintf("%s_RI_%s.png", ToolBox.folder_name, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, ...
    sprintf("%s_RI_%s.eps", ToolBox.folder_name, name)));
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
yline(vMax_mean, '--k', sprintf('%.1f mm/s', vMax_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
yline(v_mean, '--k', sprintf('%.1f mm/s', v_mean), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
yline(vMin_mean, '--k', sprintf('%.1f mm/s', vMin_mean), ...
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

box on
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])

% Export
exportgraphics(gcf, fullfile(ToolBox.path_png, sprintf("%s_PI_%s.png", ToolBox.folder_name, name)));
exportgraphics(gcf, fullfile(ToolBox.path_eps, sprintf("%s_PI_%s.eps", ToolBox.folder_name, name)));
close;

% Save image

if size(v_video, 3) > 1 % if given a video, output the image of RI / PI

    v_video_interp = v_video(:, :, sysIdxList(1):sysIdxList(2)); % not really interp but ok for now
    v_video_interp = imresize3(v_video_interp, [size(v_video_interp, 1), size(v_video_interp, 2), Ninterp]);

    vMax = mean(v_video_interp(:, :, vMax_frames_indxs), 3);
    vMin = mean(v_video_interp(:, :, vMin_frames_indxs), 3);

    % Clip RI to [0, 1] and handle NaNs
    RI = (vMax - vMin) ./ vMax;
    RI(RI > 1) = 1;
    RI(RI < 0) = 0;
    RI(isnan(RI)) = 0;

    % Clip PI to [0, -] and handle NaNs
    PI = (vMax - vMin) ./ v_mean;
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
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_RI_map_%s.png", ToolBox.folder_name, name)));

    % Display and save the PI image
    f = figure("Visible", "off");
    imagesc(PI), axis image; axis off;
    colorbar, colormap(cmap)
    title(sprintf('PI %s = %0.2f', name, PI_mean));
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_PI_map_%s.png", ToolBox.folder_name, name)));

    % Close figures
    close(f), close(fig);

else

end

if contains(name, 'artery')
    VesselName = 'arterial';
else
    VesselName = 'venous';
end

% Save json

if contains(name, 'velocity')
    ToolBox.outputs.velocity.(sprintf('mean_%s', VesselName)) = round(v_mean, 2);
    ToolBox.outputs.velocity.(sprintf('mean_%s_se', VesselName)) = round(std(signal), 2);
    ToolBox.outputs.velocity.(sprintf('systolic_%s', VesselName)) = round(vMax_mean, 2);
    ToolBox.outputs.velocity.(sprintf('systolic_%s_se', VesselName)) = round(vMax_std, 2);
    ToolBox.outputs.velocity.(sprintf('diastolic_%s', VesselName)) = round(vMin_mean, 2);
    ToolBox.outputs.velocity.(sprintf('diastolic_%s_se', VesselName)) = round(vMin_std, 2);

    % New
    if contains(name, 'vein')
        ToolBox.Outputs.add('VenousMeanVelocity', v_mean, 'mm/s', std(signal));
        ToolBox.Outputs.add('VenousMaximumVelocity', vMax_mean, 'mm/s', vMin_std);
        ToolBox.Outputs.add('VenousMinimumVelocity', vMin_mean, 'mm/s', vMax_std);
    end

end

ToolBox.outputs.indices.(sprintf('%s_RI', name)) = round(RI_mean, 2);
ToolBox.outputs.indices.(sprintf('%s_RI_se', name)) = round(dRI, 2);
ToolBox.outputs.indices.(sprintf('%s_PI', name)) = round(PI_mean, 2);
ToolBox.outputs.indices.(sprintf('%s_PI_se', name)) = round(dPI, 2);
ToolBox.outputs.indices.(sprintf('%s_PR', name)) = round(PR_mean, 2);
ToolBox.outputs.indices.(sprintf('%s_PR_se', name)) = round(dPR, 2);

% New
if contains(name, 'vein')
    ToolBox.Outputs.add('VenousResistivityIndexVelocity', RI_mean, '', dRI);
    ToolBox.Outputs.add('VenousPulsatilityIndexVelocity', PI_mean, '', dPI);
    ToolBox.Outputs.add('VenousMaxMinRatioVelocity', (vMax_mean / vMin_mean), '');
end

end
