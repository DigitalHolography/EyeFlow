function [v_RMS_video] = pulseAnalysis(f_RMS_video, maskArtery, maskVein, maskNeighbors, xy_barycenter)
% pulseAnalysis.m computes the velocities
% Inputs:
%       VIDEOS:
%   f_RMS_video     Size: numX x numY x numFrames double
%   M0_disp_video   Size: numX x numY x numFrames double
%       IMAGES:
%   f_AVG_image     Size: numX x numY double
%   maskArtery      Size: numX x numY logical
%   maskBackground  Size: numX x numY logical
%   maskVein        Size: numX x numY logical
%       TRIVIA:
%   sysIdxList:     Size: numSystoles
%
% Output:
%   v_RMS_video     Size: numX x numY x numFrames double

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;

[numX, numY, numFrames] = size(f_RMS_video);
x_c = xy_barycenter(1) / numX;
y_c = xy_barycenter(2) / numY;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
maskSection = diskMask(numX, numY, r1, r2, center = [x_c, y_c]);

maskArterySection = maskArtery & maskSection;
maskVeinSection = maskVein & maskSection;
maskVesselSection = (maskVein | maskArtery) & maskSection;

folder = 'bloodFlowVelocity';

strXlabel = 'Time(s)'; %createXlabelTime(1);
strYlabel = 'frequency (kHz)';
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
cBlack = [0 0 0];
cArtery = [255 22 18] / 255;
cVein = [18 23 255] / 255;

%% 1) Local BKG Artery and Veins %~1min

tic

diaphragmRadius = params.json.Mask.DiaphragmRadius;
maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
f_RMS_background = zeros(numX, numY, numFrames, 'single');

if veinsAnalysis
    maskVessel = maskArtery | maskVein;
else
    maskVessel = maskArtery;
end

w = params.json.PulseAnalysis.LocalBackgroundWidth;
k = params.json.Preprocess.InterpolationFactor;
bkg_scaler = params.json.PulseAnalysis.bkgScaler;

if params.json.Mask.AllNonVesselsAsBackground
    SE = strel('disk', params.json.PulseAnalysis.LocalBackgroundWidth);
    maskNeighbors = imerode(maskNeighbors, SE);

    parfor frameIdx = 1:numFrames
        f_RMS_background(:, :, frameIdx) = single(regionfill(f_RMS_video(:, :, frameIdx), ~maskNeighbors & maskDiaphragm));
    end

else

    parfor frameIdx = 1:numFrames
        f_RMS_background(:, :, frameIdx) = single(maskedAverage(f_RMS_video(:, :, frameIdx), bkg_scaler * w * 2 ^ k, maskNeighbors, maskVessel));
    end

end

imwrite(rescale(squeeze(mean(f_RMS_background, 3))), fullfile(ToolBox.path_png, folder, sprintf("%s_frequency_RMS_bkg.png", ToolBox.main_foldername)));

f_RMS_Artery = squeeze(sum(f_RMS_video .* maskArterySection, [1, 2]) / nnz(maskArterySection));
f_bkg_Artery = squeeze(sum(f_RMS_background .* maskArterySection, [1, 2]) / nnz(maskArterySection));

graphSignal('Arteries_fRMS', folder, ...
    t, f_RMS_Artery, '-', cArtery, ...
    t, f_bkg_Artery, '--', cBlack, ...
    Title = 'Average f_{RMS} in Arteries', xlabel = strXlabel, ylabel = strYlabel, ...
    Legend = {'Arteries', 'Local Background'});

fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_advanced_outputs', '.txt')), 'a');
fprintf(fileID, 'Mean fRMS difference artery : %f (kHz) \r\n', mean(f_RMS_Artery) - mean(f_bkg_Artery));
fclose(fileID);

if veinsAnalysis
    f_RMS_Vein = squeeze(sum(f_RMS_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection));
    f_bkg_Vein = squeeze(sum(f_RMS_background .* maskVeinSection, [1, 2]) / nnz(maskVeinSection));
    f_bkg_Vessel = squeeze(sum(f_RMS_background .* maskVesselSection, [1, 2]) / nnz(maskVesselSection));

    graphSignal('Veins_fRMS', folder, ...
        t, f_RMS_Vein, '-', cVein, ...
        t, f_bkg_Vein, '--', cBlack, ...
        Title = 'Average f_{RMS} in Veins', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Veins', 'Local Background'});
    fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_advanced_outputs', '.txt')), 'a');
    fprintf(fileID, 'Mean fRMS difference vein : %f (kHz) \r\n', mean(f_RMS_Vein) - mean(f_bkg_Vein));
    fclose(fileID);

    graphSignal('Vascular_fRMS', folder, ...
        t, f_RMS_Artery, '-', cArtery, ...
        t, f_RMS_Vein, '-', cVein, ...
        t, f_bkg_Vessel, '--', cBlack, ...
        Title = 'Average f_{RMS} in Vessels', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Arteries', 'Veins', 'Local Background'});
end

fprintf("    1. Local BKG Artery and Veins calculation took %ds\n", round(toc))

%% 2) Difference calculation

tic

if params.json.PulseAnalysis.DifferenceMethods == 0 %SIGNED DIFFERENCE FIRST
    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    delta_f_RMS = sign(tmp) .* sqrt(abs(tmp));
    clear tmp
elseif params.json.PulseAnalysis.DifferenceMethods == 1 % DIFFERENCE FIRST
    tmp = f_RMS_video .^ 2 - f_RMS_background .^ 2;
    tmp = tmp .* (tmp > 0);
    delta_f_RMS = sqrt(tmp);
    clear tmp
else % DIFFERENCE LAST
    delta_f_RMS = f_RMS_video - f_RMS_background;
end

% Delta f_rms plots
delta_f_RMS_Artery = delta_f_RMS .* maskArterySection;
delta_f_RMS_Artery(~maskArterySection) = NaN;
delta_f_RMS_Artery_Signal = squeeze(sum(delta_f_RMS_Artery, [1, 2], 'omitnan') / nnz(maskArterySection))';
delta_f_RMS_Artery_std = squeeze(std(delta_f_RMS_Artery, [], [1, 2], 'omitnan'))';

% Create figure for delta f_RMS in arteries
fig1 = figure;
graphSignalStd(fig1, delta_f_RMS_Artery_Signal, delta_f_RMS_Artery_std, numFrames, ...
    'frequency (kHz)', strXlabel, 'Average estimated frequency in Arteries', 'kHz');
exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_frequency_artery.png')))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_frequency_artery.eps')))

if veinsAnalysis
    delta_f_RMS_Vein = delta_f_RMS .* maskVeinSection;
    delta_f_RMS_Vein(~maskVeinSection) = NaN;
    delta_f_RMS_Vein_Signal = squeeze(sum(delta_f_RMS_Vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';
    delta_f_RMS_Vein_std = squeeze(std(delta_f_RMS_Vein, [], [1, 2], 'omitnan'))';

    fig2 = figure;
    graphSignalStd(fig2, delta_f_RMS_Vein_Signal, delta_f_RMS_Vein_std, numFrames, ...
        'frequency (kHz)', strXlabel, 'Average estimated frequency in Veins', 'kHz');
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_frequency_vein.png')))
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_frequency_vein.eps')))

    % Combined plot for arteries and veins
    graphSignal('2_Vessels_frequency', folder, ...
        t, delta_f_RMS_Artery_Signal, '-', cArtery, ...
        t, delta_f_RMS_Vein_Signal, '-', cVein, ...
        Title = 'Average estimated frequency in Arteries and Veins', xlabel = strXlabel, ylabel = 'frequency (kHz)');
else
    graphSignal('2_Arteries_frequency', folder, ...
        t, delta_f_RMS_Artery_Signal, '-', cArtery, ...
        Title = 'Average estimated frequency in Arteries', xlabel = strXlabel, ylabel = 'frequency (kHz)');
end

% Velocity calculation
scalingFactor = 1000 * 1000 * 2 * params.json.PulseAnalysis.Lambda / sin(params.json.PulseAnalysis.Phi);
v_RMS_video = scalingFactor * delta_f_RMS;

% Velocity plots
v_Artery = v_RMS_video .* maskArterySection;
v_Artery(~maskArterySection) = NaN;
v_Artery_Signal = squeeze(sum(v_Artery, [1, 2], 'omitnan') / nnz(maskArterySection))';
v_Artery_std = squeeze(std(v_Artery, [], [1, 2], 'omitnan'))';

% Create figure for velocity in arteries
fig3 = figure;
graphSignalStd(fig3, v_Artery_Signal, v_Artery_std, numFrames, ...
    'Velocity (mm/s)', strXlabel, 'Average estimated velocity in Arteries', 'mm/s');
exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_velocity_artery.png')))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_velocity_artery.eps')))

if veinsAnalysis
    v_Vein = v_RMS_video .* maskVeinSection;
    v_Vein(~maskVeinSection) = NaN;
    v_Vein_Signal = squeeze(sum(v_Vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';
    v_Vein_std = squeeze(std(v_Vein, [], [1, 2], 'omitnan'))';

    fig4 = figure;
    graphSignalStd(fig4, v_Vein_Signal, v_Vein_std, numFrames, ...
        'Velocity (mm/s)', strXlabel, 'Average estimated velocity in Veins', 'mm/s');
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_velocity_vein.png')))
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '2_velocity_vein.eps')))

    % Combined plot for arteries and veins
    graphSignal('2_Vessels_velocity', folder, ...
        t, v_Artery_Signal, '-', cArtery, ...
        t, v_Vein_Signal, '-', cVein, ...
        Title = 'Average estimated velocity in Arteries and Veins', xlabel = strXlabel, ylabel = 'Velocity (mm/s)');
else
    graphSignal('2_Arteries_velocity', folder, ...
        t, v_Artery_Signal, '-', cArtery, ...
        Title = 'Average estimated velocity in Arteries', xlabel = strXlabel, ylabel = 'Velocity (mm/s)');
end

fprintf("    2. Difference calculation took %ds\n", round(toc))

ArterialResistivityIndex(t, v_RMS_video, maskArtery .* maskSection, 'velocityArtery', folder);
ArterialResistivityIndex(t, v_RMS_video, maskVein .* maskSection, 'velocityVein', folder);

%% 3) Plots of f_RMS mean Local Background in vessels and Delta frequency in vessels and their colorbars
tic

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];

LocalBackground_in_vessels = mean(f_RMS_background, 3) .* maskVessel + ones(numX, numY) * mean(sum(f_RMS_background .* maskVessel, [1, 2]) / nnz(maskVessel), 3) .* ~maskVessel;
imagesc(LocalBackground_in_vessels);
colormap gray
title('Local Background in vessels');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS Doppler frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;

imwrite(rescale(LocalBackground_in_vessels), fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_LocalBackground_in_vessels.png')))

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
LocalBackground_colorbar = colorbar('north');
clim(range)
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
LocalBackground_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(LocalBackground_colorbar, 'Title');
titleString = 'Local Background RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_ColorBarLocalBackground_in_vessels.png')))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_ColorBarLocalBackground_in_vessels.eps')))

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];
in_vessels = mean(delta_f_RMS, 3) .* maskVessel;
imagesc(in_vessels);
colormap gray
title('Delta f in vessels');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'Delta Doppler RMS frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;
imwrite(rescale(in_vessels), fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_Df_in_vessels.png')))

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
Df_colorbar = colorbar('north');
clim(range);
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
Df_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(Df_colorbar, 'Title');
titleString = 'Delta Doppler RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_ColorBarDf_in_vessels.png')))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_ColorBarDf_in_vessels.eps')))

figure("Visible", "off")
imagesc(squeeze(mean(f_RMS_video, 3)));
colormap gray
title('RMS frequency map RAW');
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = 'RMS frequency (kHz)';
c.Label.FontSize = 12;
axis off
axis image
range(1:2) = clim;
imwrite(rescale(squeeze(mean(f_RMS_video, 3))), fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_frequency_RMS.png')), 'png');

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
f_RMS_colorbar = colorbar('north');
clim(range);
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
f_RMS_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(f_RMS_colorbar, 'Title');
titleString = 'RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_ColorbarRMSFrequency.png')))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s", ToolBox.main_foldername, '3_ColorbarRMSFrequency.eps')))

fprintf("    3. Plotting heatmaps took %ds\n", round(toc))

if exportVideos
    f_RMS_video_rescale = rescale(f_RMS_video);
    f_RMS_background_rescale = rescale(f_RMS_background);

    writeGifOnDisc(imresize(f_RMS_background_rescale, 0.5), "f_RMS_bkg")
    writeGifOnDisc(imresize(f_RMS_video_rescale, 0.5), "f_RMS")
end

clear LocalBackground_in_vessels f_RMS_background
close all

end
