function crossSectionsFigures(Q_results, name, M0_ff_video, xy_barycenter, systolesIndexes, sysIdx, diasIdx, v_video_RGB, v_mean_RGB)

% 0. Initialise Variables

ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
path_eps = ToolBox.path_eps;
params = ToolBox.getParams;
initial = name(1);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));
[numX, numY, numFrames] = size(M0_ff_video);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);

t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

v_cell = Q_results.v_cell;
v_profiles_cell = Q_results.v_profiles_cell;
dv_profiles_cell = Q_results.dv_profiles_cell;
D_cell = Q_results.D_cell;
dD_cell = Q_results.dD_cell;
locsLabel = Q_results.locsLabel;
maskLabel = Q_results.maskLabel;
A_cell = Q_results.A_cell;
Q_cell = Q_results.Q_cell;
radiusQ = Q_results.radiusQ;
radiusQSE = Q_results.radiusQSE;
branchQ = Q_results.branchQ;
labeledVessels = Q_results.labeledVessels .* Q_results.labeledVessels ~= 0;

% 1. SBlood Volume Rate Figures
tic

if ~isempty(systolesIndexes)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end

[Q_t, dQ_t] = plotRadius(radiusQ, radiusQSE, t, index_start, index_end, name);

if contains(name, 'Artery')
    ToolBox.Signals.add('ArterialVolumeRate', Q_t, 'µL/min', t, 's', dQ_t);
else
    ToolBox.Signals.add('VenousVolumeRate', Q_t, 'µL/min', t, 's', dQ_t);
end

r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
maskSection = diskMask(numX, numY, r1, r2, center = [x_c / numX y_c / numY]);
s = regionprops(labeledVessels & maskSection, 'centroid');
centroids = cat(1, s.Centroid);

graphCombined(M0_ff_video, v_video_RGB, v_mean_RGB, ...
    labeledVessels .* maskSection, ...
    Q_t, dQ_t, xy_barycenter, sprintf('%s_vr', name), ...
    'etiquettes_locs', centroids, ...
    'etiquettes_values', branchQ);

if params.json.CrossSectionsFigures.BloodFlowProfiles
    interpolatedBloodVelocityProfile(v_profiles_cell, dv_profiles_cell, sysIdx, diasIdx, name)
end

fprintf("    1. Blood Volume Rate Figures (%s) took %ds\n", name, round(toc))

% 2. Sections Figures

tic

if params.json.CrossSectionsFigures.sectionImage
    sectionImage(M0_ff_img, maskLabel, initial)
end

if params.json.CrossSectionsFigures.circleImages

    if ~isfolder(fullfile(ToolBox.path_png, 'local', 'vesselSegmentImages'))
        mkdir(fullfile(path_png, 'local'), 'vesselSegmentImages')
        mkdir(fullfile(path_eps, 'local'), 'vesselSegmentImages')
        mkdir(fullfile(path_png, 'local', 'vesselSegmentImages'), 'lumenDiameter')
        mkdir(fullfile(path_eps, 'local', 'vesselSegmentImages'), 'lumenDiameter')
        mkdir(fullfile(path_png, 'local', 'vesselSegmentImages'), 'vesselSegmentId')
        mkdir(fullfile(path_eps, 'local', 'vesselSegmentImages'), 'vesselSegmentId')
        mkdir(fullfile(path_png, 'local', 'vesselSegmentImages'), 'bloodVolumeRate')
        mkdir(fullfile(path_eps, 'local', 'vesselSegmentImages'), 'bloodVolumeRate')
        mkdir(fullfile(path_png, 'local', 'vesselSegmentImages'), 'velocity')
        mkdir(fullfile(path_eps, 'local', 'vesselSegmentImages'), 'velocity')
    end

    circleImages(M0_ff_img, xy_barycenter, A_cell, Q_cell, v_cell, maskLabel, locsLabel, name)
end

if params.json.CrossSectionsFigures.widthHistogram
    [~, avgWidth, stdWidth] = widthHistogram(D_cell, dD_cell, A_cell, name);
end

fprintf("    2. Sections Images Figures (%s) took %ds\n", name, round(toc))

% 3. Arterial Indicators

tic

if params.json.CrossSectionsFigures.strokeAndTotalVolume && ~isempty(systolesIndexes)
    strokeAndTotalVolume(Q_t, dQ_t, systolesIndexes, t, 1000, name);
end

if params.json.CrossSectionsFigures.ARIBVR
    ArterialResistivityIndex(Q_t, systolesIndexes, sprintf('BVR%s', name), 'local', dQ_t);
end

% Define parameters with uncertainties
deltaP = [6786, 280];          % Pa
avg_r = avgWidth * 1e-6 / 2;
std_r = stdWidth * 1e-6 / 2;
r = [avg_r, std_r];            % m
L = [5e-3, 1e-3];              % m

calculateHemodynamicParameters(Q_t, dQ_t, deltaP, r, L, index_start, index_end);

fprintf("    3. Arterial Indicators Images Generation (%s) took %ds\n", name, round(toc))

close all

end
