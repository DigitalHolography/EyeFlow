function [] = ArterialResistivityIndex(t, v_video, maskArtery, name, folder)

ToolBox = getGlobalToolBox;
% Color Maps
cArtery = [255 22 18] / 255;
v_video = v_video .* maskArtery;

arterial_signal = squeeze(sum(v_video .* maskArtery, [1 2])) / nnz(maskArtery)';
[~, ~, ~, ~, sysindexes, diasindexes] = compute_diasys(v_video, maskArtery);
vSys = mean(v_video(:, :, sysindexes), 3);
vDias = mean(v_video(:, :, diasindexes), 3);

v_mean = mean(arterial_signal);
vSys_mean = mean(arterial_signal(sysindexes));
vDias_mean = mean(arterial_signal(diasindexes));

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
    t, arterial_signal, '-', cArtery, ...
    ylabel = 'Velocity (mm/s)', xlabel = 'Time (s)', ...
    Title = sprintf('ARI %s = %0.2f', name, ARI_mean), ...
    yLines = [vDias_mean, vSys_mean], yLineLabels = {sprintf("%.0f mm/s", vDias_mean), sprintf("%.0f mm/s", vSys_mean)});

% API Graph
graphSignal(sprintf('API_%s', name), folder, ...
    t, arterial_signal, '-', cArtery, ...
    ylabel = 'Velocity (mm/s)', xlabel = 'Time (s)', ...
    Title = sprintf('API %s = %0.2f', name, API_mean), ...
    yLines = [vDias_mean, v_mean, vSys_mean], yLineLabels = {sprintf("%.0f mm/s", vDias_mean), sprintf("%.0f mm/s", v_mean), sprintf("%.0f mm/s", vSys_mean)});

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
end

end
