function [v_RMS_video, sysIdxList, sysIdx, diasIdx] = pulseAnalysis(f_video, maskArtery, maskVein, maskNeighbors, xy_barycenter)
% pulseAnalysis.m computes the velocities
% Inputs:
%       VIDEOS:
%   f_video     Size: numX x numY x numFrames double
%   M0_disp_video   Size: numX x numY x numFrames double
%       IMAGES:
%   f_AVG_image     Size: numX x numY double
%   maskArtery      Size: numX x numY logical
%   maskBackground  Size: numX x numY logical
%   maskVein        Size: numX x numY logical
%
% Output:
%   v_RMS_video     Size: numX x numY x numFrames double

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
veinsAnalysis = params.veins_analysis;
exportVideos = params.exportVideos;
scalingFactor = 1000 * 1000 * 2 * params.json.PulseAnalysis.Lambda / sin(params.json.PulseAnalysis.Phi);

[numX, numY, numFrames] = size(f_video);
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

% 1) Local BKG Artery and Veins %~1min

% Skeletonize and label the individual branches
skel = bwskel(maskArtery);
skel = skel & maskSection; % takes out arteries near the center
[label, n] = bwlabel(skel & ~imdilate(bwmorph(skel, 'branchpoints'), strel('disk', 2))); % labeling the skeleton with branch points off to get individual branches

% Get the label mask back to initial size
masks = zeros(numX, numY, n);

for i = 1:n % for each individual branch
    sk_mask = label == i;
    label(bwareafilt(imdilate(sk_mask, strel('disk', floor(numX * r1 / 10))) & maskArtery, 1)) = i;
    masks(:, :, i) = label == i;
end

tic

diaphragmRadius = params.json.Mask.DiaphragmRadius;
maskDiaphragm = diskMask(numX, numY, diaphragmRadius);
f_bkg = zeros(numX, numY, numFrames, 'single');

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
        f_bkg(:, :, frameIdx) = single(regionfill(f_video(:, :, frameIdx), ~maskNeighbors & maskDiaphragm));
    end

else

    parfor frameIdx = 1:numFrames
        f_bkg(:, :, frameIdx) = single(maskedAverage(f_video(:, :, frameIdx), bkg_scaler * w * 2 ^ k, maskNeighbors, maskVessel));
    end

end

frequencies_artery = zeros(numFrames, n);

for i = 1:n
    frequencies_artery(:, i) = scalingFactor * squeeze(sum((f_video - f_bkg) .* masks(:, :, i), [1, 2]) / nnz(masks(:, :, i)));
    % frequencies_artery(:, i) = wdenoise(frequencies_artery(:, i), 4);
end

imwrite(rescale(corrcoef(frequencies_artery)), fullfile(ToolBox.path_png, folder, sprintf("%s_corr_coeff_map.png", ToolBox.main_foldername)))

f_artery = squeeze(sum(f_video .* maskArterySection, [1, 2]) / nnz(maskArterySection));
f_artery_bkg = squeeze(sum(f_bkg .* maskArterySection, [1, 2]) / nnz(maskArterySection));

graphSignal('f_artery', folder, ...
    t, f_artery, '-', cArtery, ...
    t, f_artery_bkg, '--', cBlack, ...
    Title = 'Average f_{RMS} in Arteries', xlabel = strXlabel, ylabel = strYlabel, ...
    Legend = {'Arteries', 'Local Background'});

fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_advanced_outputs', '.txt')), 'a');
fprintf(fileID, 'Mean fRMS difference artery : %f (kHz) \r\n', mean(f_artery) - mean(f_artery_bkg));
fclose(fileID);

if veinsAnalysis
    f_vein = squeeze(sum(f_video .* maskVeinSection, [1, 2]) / nnz(maskVeinSection));
    f_vein_bkg = squeeze(sum(f_bkg .* maskVeinSection, [1, 2]) / nnz(maskVeinSection));
    f_vessel_bkg = squeeze(sum(f_bkg .* maskVesselSection, [1, 2]) / nnz(maskVesselSection));

    graphSignal('f_vein', folder, ...
        t, f_vein, '-', cVein, ...
        t, f_vein_bkg, '--', cBlack, ...
        Title = 'Average f_{RMS} in Veins', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Veins', 'Local Background'});
    fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'EF_advanced_outputs', '.txt')), 'a');
    fprintf(fileID, 'Mean fRMS difference vein : %f (kHz) \r\n', mean(f_vein) - mean(f_vein_bkg));
    fclose(fileID);

    graphSignal('f_vascular', folder, ...
        t, f_artery, '-', cArtery, ...
        t, f_vein, '-', cVein, ...
        t, f_vessel_bkg, '--', cBlack, ...
        Title = 'Average f_{RMS} in Vessels', xlabel = strXlabel, ylabel = strYlabel, ...
        Legend = {'Arteries', 'Veins', 'Local Background'});
end

fprintf("    1. Local BKG Artery and Veins calculation took %ds\n", round(toc))

% 2) Difference calculation

tic

if params.json.PulseAnalysis.DifferenceMethods == 0 %SIGNED DIFFERENCE FIRST
    tmp = f_video .^ 2 - f_bkg .^ 2;
    df = sign(tmp) .* sqrt(abs(tmp));
    clear tmp
elseif params.json.PulseAnalysis.DifferenceMethods == 1 % DIFFERENCE FIRST
    tmp = f_video .^ 2 - f_bkg .^ 2;
    tmp = tmp .* (tmp > 0);
    df = sqrt(tmp);
    clear tmp
else % DIFFERENCE LAST
    df = f_video - f_bkg;
end

% Outlier cleaning

df_cleaned = df; % Initialize with the original video
window_size = params.json.Preprocess.OutlierAnalysisWindow; % Compute the average profile over time (mean intensity per frame)
frame_means = squeeze(mean(df, [1 2])); % 1D array: numFrames x 1
outlier_frames_mask = isoutlier(frame_means, 'movmedian', window_size);

% Find the indices of outlier frames
outlier_indices = find(outlier_frames_mask)';

% For each outlier frame, interpolate linearly using neighboring frames
for idx = outlier_indices
    % Find the previous and next non-outlier frames
    prev_frame = find(~outlier_frames_mask(1:idx - 1), 1, 'last'); % Last non-outlier before idx
    next_frame = find(~outlier_frames_mask(idx + 1:end), 1, 'first') + idx; % First non-outlier after idx

    % Handle edge cases (e.g., first or last frame is an outlier)
    if isempty(prev_frame)
        prev_frame = next_frame; % If no previous frame, use the next frame
    end

    if isempty(next_frame)
        next_frame = prev_frame; % If no next frame, use the previous frame
    end

    if prev_frame == next_frame
        df_cleaned(:, :, idx) = df_cleaned(:, :, prev_frame);
    else
        % Linearly interpolate the outlier frame
        alpha = (idx - prev_frame) / (next_frame - prev_frame); % Interpolation weight
        df_cleaned(:, :, idx) = (1 - alpha) * df(:, :, prev_frame) + alpha * df(:, :, next_frame);
    end

end

% Delta f plots
df_artery = df .* maskArterySection;
df_artery(~maskArterySection) = NaN;
df_artery_signal = squeeze(sum(df_artery, [1, 2], 'omitnan') / nnz(maskArterySection))';
df_artery_ste = squeeze(std(df_artery, [], [1, 2], 'omitnan'))';

% Create figure for delta f in arteries
fig1 = figure("Visible", "off");
graphSignalStd(fig1, df_artery_signal, df_artery_ste, numFrames, ...
    'frequency (kHz)', strXlabel, ...
    'Average frequency in Arteries', 'kHz', 'ToolBox', ToolBox);
exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_f_artery.png", ToolBox.main_foldername)))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_f_artery.eps", ToolBox.main_foldername)))

if veinsAnalysis
    df_vein = df .* maskVeinSection;
    df_vein(~maskVeinSection) = NaN;
    df_vein_signal = squeeze(sum(df_vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';
    df_vein_std = squeeze(std(df_vein, [], [1, 2], 'omitnan'))';

    fig2 = figure("Visible", "off");
    graphSignalStd(fig2, df_vein_signal, df_vein_std, numFrames, ...
        'frequency (kHz)', strXlabel, ...
        'Average frequency in Veins', 'kHz', 'ToolBox', ToolBox);
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_f_vein.png", ToolBox.main_foldername)))
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_f_vein.eps", ToolBox.main_foldername)))

    % Combined plot for arteries and veins
    graphSignal('2_Vessels_frequency', folder, ...
        t, df_artery_signal, '-', cArtery, ...
        t, df_vein_signal, '-', cVein, ...
        Title = 'Average frequency in Arteries and Veins', xlabel = strXlabel, ylabel = 'frequency (kHz)');
else
    graphSignal('2_Arteries_frequency', folder, ...
        t, df_artery_signal, '-', cArtery, ...
        Title = 'Average frequency in Arteries', xlabel = strXlabel, ylabel = 'frequency (kHz)');
end

% Velocity calculation
v_RMS_video = scalingFactor * df;

% Velocity plots
v_artery = v_RMS_video .* maskArterySection;
v_artery(~maskArterySection) = NaN;
v_artery_signal = squeeze(sum(v_artery, [1, 2], 'omitnan') / nnz(maskArterySection))';
v_artery_ste = squeeze(std(v_artery, [], [1, 2], 'omitnan'))';

% Create figure for velocity in arteries
fig3 = figure("Visible", "off");
graphSignalStd(fig3, v_artery_signal, v_artery_ste, numFrames, ...
    'Velocity (mm/s)', strXlabel, ...
    'Average velocity in Arteries', 'mm/s', 'ToolBox', ToolBox);
ToolBox.Signals.add('ArterialVelocity', v_artery_signal, 'mm/s', t, 's', v_artery_ste);
exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_v_artery.png", ToolBox.main_foldername)))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_v_artery.eps", ToolBox.main_foldername)))

if veinsAnalysis
    v_vein = v_RMS_video .* maskVeinSection;
    v_vein(~maskVeinSection) = NaN;
    v_vein_signal = squeeze(sum(v_vein, [1, 2], 'omitnan') / nnz(maskVeinSection))';
    v_vein_ste = squeeze(std(v_vein, [], [1, 2], 'omitnan'))';

    fig4 = figure("Visible", "off");
    graphSignalStd(fig4, v_vein_signal, v_vein_ste, numFrames, ...
        'Velocity (mm/s)', strXlabel, ...
        'Average velocity in Veins', 'mm/s', 'ToolBox', ToolBox);
    ToolBox.Signals.add('VenousVelocity', v_vein_signal, 'mm/s', t, 's', v_vein_ste);
    exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_v_vein.png", ToolBox.main_foldername')))
    exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_v_vein.eps", ToolBox.main_foldername)))

    % Combined plot for arteries and veins
    graphSignal('2_Vessels_velocity', folder, ...
        t, v_artery_signal, '-', cArtery, ...
        t, v_vein_signal, '-', cVein, ...
        Title = 'Average velocity in Arteries and Veins', xlabel = strXlabel, ylabel = 'Velocity (mm/s)');
else
    graphSignal('2_Arteries_velocity', folder, ...
        t, v_artery_signal, '-', cArtery, ...
        Title = 'Average velocity in Arteries', xlabel = strXlabel, ylabel = 'Velocity (mm/s)');
end

fprintf("    2. Difference calculation took %ds\n", round(toc))

findSystoleTimer = tic;

[sysIdxList, fullPulse, sysMaxList, sysMinList] = find_systole_index(v_RMS_video, maskArtery);
[~, ~, ~, ~, sysIdx, diasIdx] = compute_diasys(v_RMS_video, maskArtery, 'bloodFlowVelocity');

% Check if the output vectors are long enough
if numel(sysIdxList) < 2 || numel(sysMaxList) < 2 || numel(sysMinList) < 2
    warning('There isnt enough systoles.');
else
    % Log systole results

    DT = ToolBox.stride / (ToolBox.fs * 1000); % Period in seconds
    HeartBeat = mean(60 ./ (diff(sysIdxList) * DT));
    HeartBeatSTE = std(60 ./ (diff(sysIdxList) * DT));
    ToolBox.Outputs.add('HeartBeat', HeartBeat, 'bpm', HeartBeatSTE);
    ToolBox.Outputs.add('SystoleIndices', sysIdxList, '');
    ToolBox.Outputs.add('MaximumSystoleIndices', sysMaxList, '');
    ToolBox.Outputs.add('MinimumDiastoleIndices', sysMinList, '');

    ToolBox.Outputs.add('TimeToMaxIncreaseSystolic', 0, 's', 0); % ref

    TimeToPeakSystole = mean((sysMaxList - sysIdxList), "omitnan") * DT;
    TimeToPeakSystoleSTE = std((sysMaxList - sysIdxList), "omitnan") * DT;
    ToolBox.Outputs.add('TimeToPeakSystole', TimeToPeakSystole, 's', TimeToPeakSystoleSTE);

    TimeToMinimumDiastole = mean((sysMinList - sysIdxList), "omitnan") * DT;
    TimeToMinimumDiastoleSTE = std((sysMinList - sysIdxList), "omitnan") * DT;
    ToolBox.Outputs.add('TimeToMinimumDiastole', TimeToMinimumDiastole, 's', TimeToMinimumDiastoleSTE);

    TimeToPeakSystoleFromMinimumDiastole = abs(TimeToMinimumDiastole) + TimeToPeakSystole;
    TimeToPeakSystoleFromMinimumDiastoleSTE = (TimeToPeakSystoleSTE + TimeToMinimumDiastoleSTE) / 2;
    ToolBox.Outputs.add('TimeToPeakSystoleFromMinimumDiastole', TimeToPeakSystoleFromMinimumDiastole, 's', TimeToPeakSystoleFromMinimumDiastoleSTE);

    Ninterp = 1000;
    interpFullPulse = interpSignal(fullPulse, sysIdxList, Ninterp);
    pMax = max(interpFullPulse);
    pMin = min(interpFullPulse);
    pRange = pMax - pMin;

    firstIndex = find(interpFullPulse - (pMin + 0.05 * pRange) < 0, 1); % Find the first index where the signal is 5 % range wise close to the min
    TimePeakToDescent = firstIndex / Ninterp * mean(diff(sysIdxList)) * DT;
    ToolBox.Outputs.add('TimePeakToDescent', TimePeakToDescent, 's');
    ToolBox.Outputs.add('TimeToDescent', TimePeakToDescent + TimeToPeakSystole, 's');

    fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_EF_main_outputs.txt')), 'a');
    fprintf(fileID, 'Heart beat: %f (bpm) \r\n', HeartBeat);
    fprintf(fileID, 'Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysIdxList), ']'));
    fprintf(fileID, 'Number of Cycles: %d \r\n', numel(sysIdxList) - 1);
    fprintf(fileID, 'Max Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMaxList), ']'));
    fprintf(fileID, 'Min Systole Indices: %s \r\n', strcat('[', sprintf("%d,", sysMinList), ']'));
    fprintf(fileID, 'Time diastolic min to systolic max derivative (ms): %f \r\n', ...
        - TimeToMinimumDiastole);
    fprintf(fileID, 'Time diastolic min to systolic max (ms): %f \r\n', ...
        TimeToPeakSystoleFromMinimumDiastole);
    fclose(fileID);

    ToolBox.outputs.HeartBeat = HeartBeat;
    ToolBox.outputs.SystoleIndices = strcat('[', sprintf("%d,", sysIdxList), ']');
    ToolBox.outputs.NumberofCycles = numel(sysIdxList) - 1;
    ToolBox.outputs.MaxSystoleIndices = strcat('[', sprintf("%d,", sysMaxList), ']');
    ToolBox.outputs.MinSystoleIndices = strcat('[', sprintf("%d,", sysMinList), ']');
    ToolBox.outputs.TimeDiastolicmintosystolicmaxderivative =- TimeToMinimumDiastole;
    ToolBox.outputs.TimeDiastolicmintosystolicmax = TimeToPeakSystoleFromMinimumDiastole;

end

fprintf("- FindSystoleIndex took: %ds\n", round(toc(findSystoleTimer)));

ArterialResistivityIndex(t, v_RMS_video, maskArtery .* maskSection, sysIdx, diasIdx, 'velocityArtery', folder);

if veinsAnalysis
    VenousResistivityIndex(t, v_RMS_video, maskVein .* maskSection, sysIdxList, 'velocityVein', folder);
end

% 3) Plots of f mean Local Background in vessels and Delta frequency in vessels and their colorbars
tic

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];

LocalBackground_in_vessels = mean(f_bkg, 3);
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

imwrite(rescale(LocalBackground_in_vessels), fullfile(ToolBox.path_png, folder, sprintf("%s_f_bkg_map.png", ToolBox.main_foldername)))

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

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_f_bkg_colorBar.png", ToolBox.main_foldername)))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_f_bkg_colorBar.eps", ToolBox.main_foldername)))

f18 = figure("Visible", "off");
f18.Position = [1100 485 350 420];
in_vessels = mean(df, 3) .* maskVessel;
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
imwrite(rescale(in_vessels), fullfile(ToolBox.path_png, folder, sprintf("%s_df_map.png", ToolBox.main_foldername)))

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

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_df_colorBar.png", ToolBox.main_foldername)))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_df_colorBar.eps", ToolBox.main_foldername)))

figure("Visible", "off")
imagesc(squeeze(mean(f_video, 3)));
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
imwrite(rescale(squeeze(mean(f_video, 3))), fullfile(ToolBox.path_png, folder, sprintf("%s_f_map.png", ToolBox.main_foldername)), 'png');

colorfig = figure("Visible", "off");
colorfig.Units = 'normalized';
colormap(c);
colormap gray
f_colorbar = colorbar('north');
clim(range);
set(gca, 'Visible', false)
set(gca, 'LineWidth', 3);
f_colorbar.Position = [0.10 0.3 0.81 0.35];
colorfig.Position(4) = 0.1000;
fontsize(gca, 15, "points");
colorTitleHandle = get(f_colorbar, 'Title');
titleString = 'RMS frequency (kHz)';
set(colorTitleHandle, 'String', titleString);

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_f_colorBar.png", ToolBox.main_foldername)))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_f_colorBar.eps", ToolBox.main_foldername)))

fprintf("    3. Plotting heatmaps took %ds\n", round(toc))

if exportVideos
    f_video_rescale = rescale(f_video);
    f_bkg_rescale = rescale(f_bkg);

    writeGifOnDisc(imresize(f_bkg_rescale, 0.5), "f_bkg")
    writeGifOnDisc(imresize(f_video_rescale, 0.5), "f")
end

cshiftn = ArterialWaveformAnalysis(v_artery_signal, v_artery_ste, t, sysIdxList, 1000, 'v_artery', ToolBox);

if veinsAnalysis
    VenousWaveformAnalysis(v_vein_signal, v_vein_ste, t, sysIdxList, 1000, 'v_vein', ToolBox)
end

close all

end
