function [] = ArterialResistivityIndex(t, v_video, mask, name, folder)

ToolBox = getGlobalToolBox;
% Color Maps
v_video = v_video .* mask;
numFrames = length(t);
strXlabel = 'Time(s)';

v_masked = v_video .* mask;
v_masked(~mask) = NaN;
v_masked_signal = squeeze(sum(v_masked, [1, 2], 'omitnan') / nnz(mask))';
v_masked_std = squeeze(std(v_masked, [],  [1, 2], 'omitnan'))';

[~, ~, ~, ~, sysindexes, diasindexes] = compute_diasys(v_video, mask);
vSys = mean(v_video(:, :, sysindexes), 3);
vDias = mean(v_video(:, :, diasindexes), 3);

v_mean = mean(v_masked_signal);
vSys_mean = mean(v_masked_signal(sysindexes));
vDias_mean = mean(v_masked_signal(diasindexes));

ARI = (vSys - vDias) ./ vSys;
ARI(ARI > 1) = 1;
ARI(ARI < 0) = 0;
ARI(isnan(ARI)) = 0;

ARI_mean = (vSys_mean - vDias_mean) ./ vSys_mean;

API = (vSys - vDias) ./ v_mean;
API(API < 0) = 0;
API(isnan(API)) = 0;

API_mean = (vSys_mean - vDias_mean) ./ v_mean;

% ARI Graph
graphSignal(sprintf('ARI_%s', name), folder, ...
    t, v_masked_signal, '-', [0 0 0], ...
    ylabel = 'Velocity (mm/s)', xlabel = 'Time (s)', ...
    Title = sprintf('ARI %s = %0.2f', name, ARI_mean), ...
    Fontsize = 12, ...
    yLines = [vDias_mean, vSys_mean], yLineLabels = {sprintf("%.0f mm/s", vDias_mean), sprintf("%.0f mm/s", vSys_mean)});

% API Graph
graphSignal(sprintf('API_%s', name), folder, ...
    t, v_masked_signal, '-', [0 0 0], ...
    ylabel = 'Velocity (mm/s)', xlabel = 'Time (s)', ...
    Title = sprintf('API %s = %0.2f', name, API_mean), ...
    Fontsize = 12, ...
    yLines = [vDias_mean, v_mean, vSys_mean], yLineLabels = {sprintf("%.0f mm/s", vDias_mean), sprintf("%.0f mm/s", v_mean), sprintf("%.0f mm/s", vSys_mean)});

fig = figure;

graphSignalStd(fig, v_masked_signal, v_masked_std, numFrames, ...
    'Velocity (mm/s)', strXlabel, sprintf('Average estimated velocity in %s', name), 'mm/s');

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4)])

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_ARI_graph_std_%s.png", ToolBox.main_foldername, name)));

% Save image

if size(v_video, 3) > 1 % if given a video, output the image of ARI / API
    % Generate colormap
    [cmapARI] = cmapLAB(256, [0 0 0], 0, [1 0 0], 1/3, [1 1 0], 2/3, [1 1 1], 1);

    % Create RGB images for visualization

    % Display and save the ARI image
    fig = figure("Visible", "off");
    imagesc(ARI), axis image; axis off;
    colorbar, colormap(cmapARI)
    title(sprintf('ARI %s = %0.2f', name, ARI_mean));
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_ARI_%s.png", ToolBox.main_foldername, name)));

    % Display and save the API image
    f = figure("Visible", "off");
    imagesc(API), axis image; axis off;
    colorbar, colormap(cmapARI)
    title(sprintf('API %s = %0.2f', name, API_mean));
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_API_%s.png", ToolBox.main_foldername, name)));

    % Close figures
    close(f), close(fig);

else

    % Save txt
    fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_main_outputs', '.txt')), 'a');

    if strcmp(name, 'velocity')
        fprintf(fileID, 'Mean Velocity artery : %f (mm/s) \r\n', v_mean);
        fprintf(fileID, 'Max Velocity artery : %f (mm/s) \r\n', vSys_mean);
        fprintf(fileID, 'Min Velocity artery : %f (mm/s) \r\n', vDias_mean);
    end

    fprintf(fileID, 'Arterial Resistivity Index (%s) : %f  \r\n', name, ARI);
    fprintf(fileID, 'Arterial Pulsatility Index (%s) : %f  \r\n', name, API);
    fclose(fileID);

    % Save json

    if strcmp(name, 'velocity')
        ToolBox.outputs.MeanVelocityArtery = v_mean;
        ToolBox.outputs.SysVelocityArtery = vSys_mean;
        ToolBox.outputs.DiaVelocityArtery = vDias_mean;
    end

    ToolBox.outputs.(sprintf('ARI(%s)',name)) = ARI;
    ToolBox.outputs.(sprintf('API(%s)',name)) = API;
end

end
