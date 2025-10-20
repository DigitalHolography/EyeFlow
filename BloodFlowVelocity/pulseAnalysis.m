function [v_RMS_video, v_video_RGB, v_mean_RGB] = pulseAnalysis(f_video, M0_ff)
% pulseAnalysis.m computes the flow velocities from Doppler data
% Inputs:
%       VIDEOS:
%   f_video         Size: numX x numY x numFrames double (Doppler Data)
%   M0_ff     Size: numX x numY x numFrames double (Display Data)
%       IMAGES:
%   maskArtery      Size: numX x numY logical
%   maskBackground  Size: numX x numY logical
%   maskVein        Size: numX x numY logical
%   xy_barycenter   Size: 1 x 2 double
%
% Output:
%   v_RMS_video     Size: numX x numY x numFrames double

% Initial Setup
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;
save_figures = params.json.save_figures;

maskArtery = ToolBox.Cache.maskArtery;
maskVein = ToolBox.Cache.maskVein;
maskNeighbors = ToolBox.Cache.maskNeighbors;
xy_barycenter = ToolBox.Cache.xy_barycenter;

% Validating inputs
if ~any(maskArtery)
    error("Given Mask Artery is empty.")
end

% Constants
% 1000 -> m to mm conversion
% 1000 -> kHz to Hz conversion

scalingFactor = 1000 * 1000 * 2 * params.json.PulseAnalysis.Lambda / sin(params.json.PulseAnalysis.Phi);
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
[numX, numY, numFrames] = size(f_video);

% Section 1: Background Calculation

tic

% Create section mask
x_c = xy_barycenter(1) / numX;
y_c = xy_barycenter(2) / numY;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
maskSection = diskMask(numX, numY, r1, r2, 'center', [x_c, y_c]);

% Create vessel masks
maskArterySection = maskArtery & maskSection;
maskVeinSection = maskVein & maskSection;
maskVesselSection = (maskVein | maskArtery) & maskSection;

% Validating inputs
if ~any(maskArterySection)
    error("Given Mask Artery has no part within the current section.")
end

if ~any(maskVesselSection)
    error("Given Mask Vein has no part within the current section.")
end

f_bkg = zeros(numX, numY, numFrames, 'single');

% Determine vessel mask based on analysis type
if veinsAnalysis
    maskVessel = maskArtery | maskVein;
else
    maskVessel = maskArtery;
end

% Background calculation parameters
w = params.json.PulseAnalysis.LocalBackgroundWidth;
k = params.json.Preprocess.InterpolationFactor;
bkg_scaler = params.json.PulseAnalysis.bkgScaler;

% Calculate background
if params.json.Mask.AllNonVesselsAsBackground
    SE = strel('disk', params.json.PulseAnalysis.LocalBackgroundWidth);

    parfor frameIdx = 1:numFrames
        f_bkg(:, :, frameIdx) = single(regionfill(f_video(:, :, frameIdx), imdilate(maskVessel, SE)));
    end

else

    parfor frameIdx = 1:numFrames
        f_bkg(:, :, frameIdx) = single(maskedAverage(f_video(:, :, frameIdx), bkg_scaler * w * 2 ^ k, maskNeighbors, maskVessel));
    end

end

% Calculate and plot artery signals
if save_figures

    % Time vector for plotting
    t = ToolBox.Cache.t;

    % Colors for plotting
    cBlack = [0 0 0];
    cArtery = [255 22 18] / 255;
    cVein = [18 23 255] / 255;

    % Artery signals
    f_artery = squeeze(sum(f_video .* maskArterySection, [1, 2]) / nnz(maskArterySection));
    f_artery_bkg = squeeze(sum(f_bkg .* maskArterySection, [1, 2]) / nnz(maskArterySection));

    graphSignal('f_artery', ...
        t, f_artery, '-', cArtery, ...
        t, f_artery_bkg, '--', cBlack, ...
        'xlabel', 'Time(s)', 'ylabel', 'frequency (kHz)', ...
        'Legend', {'arteries', 'background'});

    if veinsAnalysis
        % Vein signals
        f_vein = squeeze(sum(f_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection));
        f_vein_bkg = squeeze(sum(f_bkg .* maskVeinSection, [1, 2]) / nnz(maskVeinSection));
        f_vessel_bkg = squeeze(sum(f_bkg .* maskVesselSection, [1, 2]) / nnz(maskVesselSection));

        graphSignal('f_vein', ...
            t, f_vein, '-', cVein, ...
            t, f_vein_bkg, '--', cBlack, ...
            'xlabel', 'Time(s)', 'ylabel', 'frequency (kHz)', ...
            'Legend', {'veins', 'background'});

        graphSignal('f_vessel', ...
            t, f_artery, '-', cArtery, ...
            t, f_vein, '-', cVein, ...
            t, f_vessel_bkg, '--', cBlack, ...
            'xlabel', 'Time(s)', 'ylabel', 'frequency (kHz)', ...
            'Legend', {'arteries', 'veins', 'background'});
    end

end

fprintf("    1. Background calculation took %ds\n", round(toc));

% Section 2: Difference Calculation and Velocity Computation

tic;

% Calculate difference based on selected method
switch params.json.PulseAnalysis.DifferenceMethods
    case 0 % SIGNED DIFFERENCE FIRST
        tmp = f_video .^ 2 - f_bkg .^ 2;
        df = sign(tmp) .* sqrt(abs(tmp));
    case 1 % DIFFERENCE FIRST
        tmp = f_video .^ 2 - f_bkg .^ 2;
        tmp = tmp .* (tmp > 0);
        df = sqrt(tmp);
    otherwise % DIFFERENCE LAST
        df = f_video - f_bkg;
end

% Clean up temporary variable
if exist('tmp', 'var'), clear tmp; end

% Process artery signals and detect outliers
df_artery = df .* maskArterySection;
df_artery_signal = squeeze(sum(df_artery, [1, 2], 'omitnan') / nnz(maskArterySection))';
outlier_frames_mask = isoutlier(df_artery_signal, "movmedian", 10);

% Process vein signals if enabled
if veinsAnalysis
    df_vein = df .* maskVeinSection;
    df_vein_signal = squeeze(sum(df_vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';
    outlier_frames_mask = outlier_frames_mask | isoutlier(df_vein_signal, "movmedian", 10);
end

% Interpolate outlier frames
df = interpolateOutlierFrames(df, outlier_frames_mask');

if save_figures
    % Recalculate signals after interpolation
    df_artery = df .* maskArterySection;
    df_artery(~maskArterySection) = NaN;
    df_artery_signal = squeeze(sum(df_artery, [1, 2], 'omitnan') / nnz(maskArterySection))';

    % Plot signals
    if veinsAnalysis
        df_vein = df .* maskVeinSection;
        df_vein(~maskVeinSection) = NaN;
        df_vein_signal = squeeze(sum(df_vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';

        graphSignal('df_vessel', ...
            t, df_artery_signal, '-', cArtery, ...
            t, df_vein_signal, '-', cVein, ...
            'xlabel', 'Time(s)', 'ylabel', 'frequency (kHz)');
    else
        graphSignal('df_artery', ...
            t, df_artery_signal, '-', cArtery, ...
            'xlabel', 'Time(s)', 'ylabel', 'frequency (kHz)');
    end

end

% Calculate velocity
v_RMS_video = scalingFactor * df;

% Process artery velocity signals
v_artery = v_RMS_video .* maskArterySection;
v_artery(~maskArterySection) = NaN;
v_artery_signal = squeeze(sum(v_artery, [1, 2], 'omitnan') / nnz(maskArterySection))';

if veinsAnalysis
    % Process vein velocity signals
    v_vein = v_RMS_video .* maskVeinSection;
    v_vein(~maskVeinSection) = NaN;
    v_vein_signal = squeeze(sum(v_vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';
end

% Filter signals
[b, a] = butter(4, 15 / (fs / 2), 'low');
v_artery_signal = filtfilt(b, a, v_artery_signal);

if veinsAnalysis
    v_vein_signal = filtfilt(b, a, v_vein_signal);
end

if save_figures
    % Plot velocity signals
    if veinsAnalysis
        graphSignal('v_vessel', ...
            t, v_artery_signal, '-', cArtery, ...
            t, v_vein_signal, '-', cVein, ...
            'Title', 'average velocity in arteries and veins', 'xlabel', 'Time(s)', 'ylabel', 'Velocity (mm/s)');
        ToolBox.Output.Signals.add('ArterialVelocity', v_artery_signal, 'mm/s', t, 's');
        ToolBox.Output.Signals.add('VenousVelocity', v_vein_signal, 'mm/s', t, 's');
    else
        graphSignal('v_artery', ...
            t, v_artery_signal, '-', cArtery, ...
            'Title', 'average velocity in arteries', 'xlabel', 'Time(s)', 'ylabel', 'Velocity (mm/s)');
        ToolBox.Output.Signals.add('ArterialVelocity', v_artery_signal, 'mm/s', t, 's');
    end

end

fprintf("    2. Difference calculation and velocity computation took %ds\n", round(toc));

% Section 3: Systole/Diastole Analysis

tic;

% Find systole indices
if veinsAnalysis
    [sys_idx_list, pulse_artery, sys_max_list, sys_min_list] = find_systole_index(v_artery_signal, 'pulseVein', v_vein_signal);
else
    [sys_idx_list, pulse_artery, sys_max_list, sys_min_list] = find_systole_index(v_artery_signal);
end

[M0_Systole_img, M0_Diastole_img, sys_idx, dias_idx] = compute_diasys(v_RMS_video, maskArterySection);

if save_figures
    v_RMS_img = mean(v_RMS_video, 3, 'omitnan');
    diasys_diff = M0_Systole_img - M0_Diastole_img;
    RGBdiasys = labDuoImage(rescale(v_RMS_img), diasys_diff);
    saveMaskImage(RGBdiasys, 'vessel_20_diasys_diff.png', isStep = true);
end

% Process heart beat data if enough cycles detected
if numel(sys_idx_list) >= 2 && numel(sys_max_list) >= 2 && numel(sys_min_list) >= 2
    DT = ToolBox.stride / (ToolBox.fs * 1000); % Period in seconds
    heartRates = 60 ./ (diff(sys_idx_list) * DT);

    % Calculate statistics
    HeartBeat = mean(heartRates);
    HeartBeatSTE = std(heartRates);
    TimeToPeakSystole = mean((sys_max_list - sys_idx_list'), "omitnan") * DT;
    TimeToPeakSystoleSTE = std((sys_max_list - sys_idx_list'), "omitnan") * DT;
    TimeToMinimumDiastole = mean((sys_min_list - sys_idx_list'), "omitnan") * DT;
    TimeToMinimumDiastoleSTE = std((sys_min_list - sys_idx_list'), "omitnan") * DT;
    TimeToPeakSystoleFromMinimumDiastole = abs(TimeToMinimumDiastole) + TimeToPeakSystole;
    TimeToPeakSystoleFromMinimumDiastoleSTE = (TimeToPeakSystoleSTE + TimeToMinimumDiastoleSTE) / 2;

    % Interpolate full pulse for additional analysis
    Ninterp = 1000;
    interpFullPulse = interpSignal(pulse_artery, sys_idx_list, Ninterp);
    pMax = max(interpFullPulse);
    pMin = min(interpFullPulse);
    pRange = pMax - pMin;

    % Find descent time
    firstIndex = find(interpFullPulse - (pMin + 0.05 * pRange) < 0, 1);
    TimePeakToDescent = firstIndex / Ninterp * mean(diff(sys_idx_list')) * DT;

    % Store output
    ToolBox.Output.add('HeartBeat', HeartBeat, 'bpm', HeartBeatSTE);
    ToolBox.Output.add('SystoleIndices', sys_idx_list, '');
    ToolBox.Output.add('ArterySystoleMaxIndices', sys_max_list, '');
    ToolBox.Output.add('ArteryDiastoleMinIndices', sys_min_list, '');
    ToolBox.Output.add('ArteryTimeToMaxIncrease', 0, 's', 0);
    ToolBox.Output.add('ArteryTimeToPeakSystole', TimeToPeakSystole, 's', TimeToPeakSystoleSTE);
    ToolBox.Output.add('ArteryTimeToMinDiastole', TimeToMinimumDiastole, 's', TimeToMinimumDiastoleSTE);
    ToolBox.Output.add('ArteryTimeToPeakSystoleFromDiastole', TimeToPeakSystoleFromMinimumDiastole, 's', TimeToPeakSystoleFromMinimumDiastoleSTE);
    ToolBox.Output.add('ArteryTimeToDescent', TimePeakToDescent, 's');
    ToolBox.Output.add('ArteryTimePeakToDescent', TimePeakToDescent + TimeToPeakSystole, 's');

    % Log detailed results
    logDetailedResults(ToolBox, HeartBeat, sys_idx_list, sys_max_list, sys_min_list, ...
        TimeToMinimumDiastole, TimeToPeakSystoleFromMinimumDiastole);
else
    warning('There isn''t enough systoles for analysis.');
end

fprintf("    3. Systole/diastole analysis took %ds\n", round(toc));

% Section 4: Resistivity Index and Waveform Analysis
tic;

% Calculate arterial resistivity index
ArterialResistivityIndex(v_artery_signal, sys_idx_list, 'v_artery', ForceFigure = true);

% Vein analysis if enabled
if veinsAnalysis
    ArterialResistivityIndex(v_vein_signal, sys_idx_list, 'v_vein', ForceFigure = true);
end

% Perform waveform analysis
v_artery_interp = ArterialWaveformAnalysis(v_artery_signal, sys_idx_list, 128, 'v_artery');

if veinsAnalysis
    v_vein_interp = VenousWaveformAnalysis(v_vein_signal, sys_idx_list, 128, 'v_vein');
end

if veinsAnalysis
    % Calculate correlation between artery and vein signals
    arterial_venous_correlation(v_artery_signal, -v_vein_signal);
end

if veinsAnalysis
    % Calculate correlation between artery and vein signals
    arterial_venous_delay(v_artery_interp, v_vein_interp);
end

fprintf("    4. Resistivity and waveform analysis took %ds\n", round(toc()));

% Section 5: Visualization and Output Generation

tic;

% background in vessels
if save_figures

    LocalBackground_in_vessels = mean(f_bkg, 3);
    createHeatmap(LocalBackground_in_vessels, 'background in vessels', ...
        'background RMS frequency (kHz)', fullfile(ToolBox.path_png, sprintf("%s_f_bkg_map.png", ToolBox.folder_name)));

    % Delta f in vessels
    in_vessels = mean(df, 3) .* maskVesselSection;
    createHeatmap(in_vessels, 'Delta f in vessels', ...
        'Delta Doppler RMS frequency (kHz)', fullfile(ToolBox.path_png, sprintf("%s_df_map_vessel.png", ToolBox.folder_name)));

    % Delta f
    in_vessels = mean(df, 3) .* maskVesselSection + (mean(f_video, 3) - mean(f_video, 'all')) .* ~maskVesselSection;
    createHeatmap(in_vessels, 'Delta f in vessels', ...
        'Delta Doppler RMS frequency (kHz)', fullfile(ToolBox.path_png, sprintf("%s_df_map.png", ToolBox.folder_name)));

    velocityIm(mean(df, 3) .* maskVesselSection, maskArtery | maskVein, turbo, 'df_vessel', colorbarOn = true, LabelName = 'kHz');

    % Raw RMS frequency map
    raw_map = squeeze(mean(f_video, 3));
    createHeatmap(raw_map, 'RMS frequency map RAW', ...
        'RMS frequency (kHz)', fullfile(ToolBox.path_png, sprintf("%s_f_map.png", ToolBox.folder_name)));

end

% Export videos if enabled
if exportVideos
    f_video_rescale = rescale(f_video);
    f_bkg_rescale = rescale(f_bkg);
    writeGifOnDisc(imresize(f_bkg_rescale, 0.5), "f_bkg");
    writeGifOnDisc(imresize(f_video_rescale, 0.5), "f");
end

v_video_RGB = zeros(numX, numY, 3, numFrames, 'single');
v_mean_RGB = zeros(numX, numY, 3, 'single');

if save_figures
    % Precompute masks for visualization
    maskAV = maskArtery & maskVein;
    maskArterySection = maskArtery & maskSection & ~maskAV;
    maskVeinSection = maskVein & maskSection & ~maskAV;

    % Generate flow maps
    [v_video_RGB, v_mean_RGB] = flowMap(v_RMS_video, maskSection, maskArtery, maskVein, M0_ff, xy_barycenter, ToolBox);

    % Generate histograms
    histoVideoArtery = VelocityHistogram(v_RMS_video, maskArterySection, 'artery');

    if veinsAnalysis
        histoVideoVein = VelocityHistogram(v_RMS_video, maskVeinSection, 'vein');

        % Generate combined visualizations
        createCombinedVisualizations(v_mean_RGB, histoVideoArtery, veinsAnalysis, histoVideoVein, ...
            v_video_RGB, numFrames, exportVideos, ToolBox);
    end

end

fprintf("    5. Visualization and output generation took %ds\n", round(toc));

close all

% Save Intermediate Results in cache
ToolBox.Cache.sysIdxList = sys_idx_list;
ToolBox.Cache.sysIdx = sys_idx;
ToolBox.Cache.diasIdx = dias_idx;

end

function logDetailedResults(ToolBox, HeartBeat, sysIdxList, sysMaxList, sysMinList, ...
    TimeToMinimumDiastole, TimeToPeakSystoleFromMinimumDiastole)
% Helper function to log detailed results
fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.folder_name, '_main_outputs.txt')), 'a');
fprintf(fileID, 'Heart beat: %f (bpm) \r\n', HeartBeat);
fprintf(fileID, 'Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysIdxList), ']'));
fprintf(fileID, 'Number of Cycles: %d \r\n', numel(sysIdxList) - 1);
fprintf(fileID, 'Max Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMaxList), ']'));
fprintf(fileID, 'Min Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMinList), ']'));
fprintf(fileID, 'Time diastolic min to systolic max derivative (ms): %f \r\n', -TimeToMinimumDiastole);
fprintf(fileID, 'Time diastolic min to systolic max (ms): %f \r\n', TimeToPeakSystoleFromMinimumDiastole);
fclose(fileID);
end

function createHeatmap(data, titleStr, cbarStr, filename)
% Helper function to create a single heatmap
f = figure("Visible", "off");
f.Position = [1100 485 350 420];
imagesc(data);
colormap gray
title(titleStr);
fontsize(gca, 12, "points");
set(gca, 'LineWidth', 2);
c = colorbar('southoutside');
c.Label.String = cbarStr;
c.Label.FontSize = 12;
axis off
axis image
rangeVals = clim;

imwrite(rescale(data), filename);

% Create colorbar figure
colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap gray
cbar = colorbar('north');
clim(rangeVals);
set(gca, 'Visible', false);
set(gca, 'LineWidth', 3);
cbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(cbar, 'Title');
set(colorTitleHandle, 'String', cbarStr);

[pathstr, name, ext] = fileparts(filename);
exportgraphics(gca, fullfile(pathstr, strrep(name, '_map', '_colorBar') + ext));
exportgraphics(gca, fullfile(strrep(pathstr, 'png', 'eps'), strrep(name, '_map', '_colorBar') + ".eps"));

close([f, colorfig]);
end

function createCombinedVisualizations(v_mean_RGB, histoVideoArtery, veinsAnalysis, histoVideoVein, ...
    v_video_RGB, numFrames, exportVideos, ToolBox)
% Helper function to create combined visualizations

% Determine sizes
[numX_fig, numY_fig, ~, ~] = size(histoVideoArtery);

% Create averaged visualization
if veinsAnalysis
    v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [numX_fig * 2 numX_fig * 2 3]));
    combinedImg = cat(2, v_mean_RGB4Gif, cat(1, mat2gray(histoVideoArtery(:, :, :, end)), mat2gray(histoVideoVein(:, :, :, end))));
else
    v_mean_RGB4Gif = rescale(imresize3(v_mean_RGB, [numY_fig numY_fig 3]));
    combinedImg = cat(1, v_mean_RGB4Gif, mat2gray(histoVideoArtery(:, :, :, end)));
end

imwrite(combinedImg, fullfile(ToolBox.path_png, sprintf("%s_%s", ToolBox.folder_name, 'AVGflowVideoCombined.png')));

% Create video visualization if exporting
if exportVideos

    if veinsAnalysis
        v_video_RGB4Gif = zeros(numX_fig * 2, numX_fig * 2, 3, numFrames, 'like', v_video_RGB);
        v_video_RGB4Gif(:, :, 1, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 1, :)), [numX_fig * 2 numX_fig * 2 numFrames]));
        v_video_RGB4Gif(:, :, 2, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 2, :)), [numX_fig * 2 numX_fig * 2 numFrames]));
        v_video_RGB4Gif(:, :, 3, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 3, :)), [numX_fig * 2 numX_fig * 2 numFrames]));
        combinedGifs = cat(2, v_video_RGB4Gif, cat(1, mat2gray(histoVideoArtery), mat2gray(histoVideoVein)));
    else
        v_video_RGB4Gif = zeros(numY_fig, numY_fig, 3, numFrames, 'like', v_video_RGB);
        v_video_RGB4Gif(:, :, 1, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 1, :)), [numY_fig numY_fig numFrames]));
        v_video_RGB4Gif(:, :, 2, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 2, :)), [numY_fig numY_fig numFrames]));
        v_video_RGB4Gif(:, :, 3, :) = mat2gray(imresize3(squeeze(v_video_RGB(:, :, 3, :)), [numY_fig numY_fig numFrames]));
        combinedGifs = cat(1, v_video_RGB4Gif, mat2gray(histoVideoArtery));
    end

    writeGifOnDisc(mat2gray(combinedGifs), "velocityHistogramCombined", 0.04);
end

end
