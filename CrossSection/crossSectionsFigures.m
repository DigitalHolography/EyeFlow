function crossSectionsFigures(results, name, M0_ff, v_video_RGB, v_mean_RGB)

% 0. Initialise Variables

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
path_png = ToolBox.path_png;
path_eps = ToolBox.path_eps;

% Retrieve cached variables
xy_barycenter = ToolBox.Cache.xy_barycenter;
systolesIndexes = ToolBox.Cache.sysIdxList;
sysIdx = ToolBox.Cache.sysIdx;
diasIdx = ToolBox.Cache.diasIdx;

initial = name(1);
M0_ff = rescale(M0_ff);
M0_ff_img = rescale(mean(M0_ff, 3));
[numX, numY, numFrames] = size(M0_ff);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);

t = ToolBox.Cache.t;

A_cell = results.A_cell;
D_cell = results.D_cell;
locsLabel = results.locsLabel;
maskLabel = results.maskLabel;
Q_cell = results.Q_cell;
v_cell = results.v_cell;
v_profiles_cell = results.v_profiles_cell;
radius_Q = results.radius_Q;
branch_Q = results.branch_Q;
radius_v = results.radius_v;
branch_v = results.branch_v;

% Standard errors
D_SE_cell = results.D_SE_cell;
v_SE_profiles_cell = results.v_SE_profiles_cell;
radius_Q_SE = results.radius_Q_SE;
branch_Q_SE = results.branch_Q_SE;
radius_v_SE = results.radius_v_SE;
branch_v_SE = results.branch_v_SE;

labeledVessels = results.labeledVessels .* results.labeledVessels ~= 0;
histo_v_cell = results.histo_v_cell;

% 1. Flow Rate Figures
tic

if ~isempty(systolesIndexes)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end

[Q_t, Q_SE_t] = plotRadius(radius_Q, radius_Q_SE, t, index_start, index_end, name, 'flux');
[v_t, v_SE_t] = plotRadius(radius_v, radius_v_SE, t, index_start, index_end, name, 'velocity');

plotBranch(branch_Q, branch_Q_SE, t, index_start, index_end, name, 'flux');
plotBranch(branch_v, branch_v_SE, t, index_start, index_end, name, 'velocity');

if contains(name, 'artery')
    ToolBox.Output.Signals.add('ArterialVolumeRate', Q_t, 'µL/min', t, 's', Q_SE_t);
elseif contains(name, 'vein')
    ToolBox.Output.Signals.add('VenousVolumeRate', Q_t, 'µL/min', t, 's', Q_SE_t);
end

if contains(name, 'artery')
    ToolBox.Output.Signals.add('ArterialVelocity', v_t, 'mm/s', t, 's', v_SE_t);
elseif contains(name, 'vein')
    ToolBox.Output.Signals.add('VenousVelocity', v_t, 'mm/s', t, 's', v_SE_t);
end

r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
maskSection = diskMask(numX, numY, r1, r2, center = [x_c / numX y_c / numY]);
s = regionprops(labeledVessels & maskSection, 'centroid');
centroids = cat(1, s.Centroid);

graphCombined(M0_ff, v_video_RGB, v_mean_RGB, ...
    labeledVessels .* maskSection, ...
    Q_t, Q_SE_t, xy_barycenter, sprintf('vr_%s', name), ...
    'etiquettes_locs', centroids, ...
    'etiquettes_values', branch_Q);

fprintf("    1. Flow Rate Figures (%s) took %ds\n", name, round(toc))

tic

if params.json.CrossSectionsFigures.BloodFlowProfiles
    interpolatedBloodVelocityProfile(v_profiles_cell, v_SE_profiles_cell, sysIdx, diasIdx, name)
end

if params.json.CrossSectionsFigures.BloodFlowHistograms
    histogramPatchVelocities(histo_v_cell, name, locsLabel, mean(M0_ff, 3))
end

if params.json.CrossSectionsFigures.BloodFlowProfilesWomersleyOverlay
    profilePatchWomersley(v_profiles_cell, name, locsLabel, mean(M0_ff, 3))
end

if params.json.CrossSectionsFigures.BloodFlowProfilesOverlay
    profilePatchVelocities(v_profiles_cell, name, locsLabel, mean(M0_ff, 3))
end

fprintf("    1.(bis) optional Flow Rate Figures (interpolated velocity profiles / Histograms / Profiles Overlay) (%s) took %ds\n", name, round(toc))

% 2. Sections Figures

tic

if params.json.CrossSectionsFigures.sectionImage
    sectionImage(M0_ff_img, maskLabel, initial)
end

if params.json.CrossSectionsFigures.circleImages

    if ~isfolder(fullfile(path_png, 'vesselSegmentImages'))
        mkdir(fullfile(path_png, 'vesselSegmentImages'))
        mkdir(fullfile(path_eps, 'vesselSegmentImages'))
        mkdir(fullfile(path_png, 'vesselSegmentImages', 'lumenDiameter'))
        mkdir(fullfile(path_eps, 'vesselSegmentImages', 'lumenDiameter'))
        mkdir(fullfile(path_png, 'vesselSegmentImages', 'vesselSegmentId'))
        mkdir(fullfile(path_eps, 'vesselSegmentImages', 'vesselSegmentId'))
        mkdir(fullfile(path_png, 'vesselSegmentImages', 'bloodVolumeRate'))
        mkdir(fullfile(path_eps, 'vesselSegmentImages', 'bloodVolumeRate'))
        mkdir(fullfile(path_png, 'vesselSegmentImages', 'velocity'))
        mkdir(fullfile(path_eps, 'vesselSegmentImages', 'velocity'))
    end

    circleImages(M0_ff_img, xy_barycenter, A_cell, Q_cell, v_cell, maskLabel, locsLabel, name)
end

if params.json.CrossSectionsFigures.widthHistogram
    [D_mid, ~, D_std] = widthHistogram(D_cell, D_SE_cell, A_cell, name);
end

fprintf("    2. Sections Images Figures (%s) took %ds\n", name, round(toc))

% 3. Arterial Indicators

tic

if params.json.CrossSectionsFigures.strokeAndTotalVolume && ~isempty(systolesIndexes)
    strokeAndTotalVolume(Q_t, Q_SE_t, systolesIndexes, 1000, name);
end

if params.json.CrossSectionsFigures.ARIBVR
    ArterialResistivityIndex(Q_t, systolesIndexes, sprintf('BVR%s', name), Q_SE_t);
end

% Define parameters with uncertainties
% deltaP = [6176, 280]; % Pa
deltaP = [1000, 100]; % (ONLY IDEAL L)
avg_r = D_mid * 1e-6/2;
std_r = D_std * 1e-6/2;
r = [avg_r, std_r]; % m
L = [5e-3, 1e-3]; % m (ONLY IDEAL L)
N = size(branch_Q, 1);

if params.json.CrossSectionsFigures.hemodynamicParameters
    calculateHemodynamicParameters(Q_t, Q_SE_t, deltaP, r, L, index_start, index_end, N);
end

fprintf("    3. Vascular Indicators Images Generation (%s) took %ds\n", name, round(toc))

close all

end
