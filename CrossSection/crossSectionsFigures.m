function crossSectionsFigures(Q_results, mask, name, M0_ff_video, xy_barycenter, systolesIndexes, sysIdx, diasIdx)

% 0. Initialise Variables

ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
path_eps = ToolBox.path_eps;
params = ToolBox.getParams;
initial = name(1);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));
[~, ~, numFrames] = size(M0_ff_video);
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

numSections = Q_results.numSections;
locs = Q_results.locs;
v_cell = Q_results.v_cell;
v_profiles_cell = Q_results.v_profiles_cell;
dv_profiles_cell = Q_results.dv_profiles_cell;
D_cell = Q_results.D_cell;
dD_cell = Q_results.dD_cell;
mask_mat = Q_results.mask_mat;
area_mat = Q_results.area_mat;
Q_cell = Q_results.Q_cell;
Q_mat = Q_results.Q_mat;
dQ_mat = Q_results.dQ_mat;

% 1. Sections Image
tic

if ~isempty(systolesIndexes)
    index_start = systolesIndexes(1);
    index_end = systolesIndexes(end);
else
    index_start = 1;
    index_end = numFrames;
end

if params.json.CrossSectionsFigures.sectionImage
    sectionImage(M0_ff_img, mask_mat, initial)
end

if params.json.CrossSectionsFigures.circleImages
    
    if ~isfolder(fullfile(ToolBox.path_png, 'crossSectionsAnalysis', 'sectionsImages'))
        mkdir(fullfile(path_png, 'crossSectionsAnalysis'), 'sectionsImages')
        mkdir(fullfile(path_eps, 'crossSectionsAnalysis'), 'sectionsImages')
        mkdir(fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages'), 'widths')
        mkdir(fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages'), 'widths')
        mkdir(fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages'), 'num')
        mkdir(fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages'), 'num')
        mkdir(fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages'), 'bvr')
        mkdir(fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages'), 'bvr')
        mkdir(fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages'), 'vel')
        mkdir(fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages'), 'vel')
    end
    
    circleImages(M0_ff_img, xy_barycenter, area_mat, Q_cell, v_cell, mask_mat, locs, name)
end

if params.json.CrossSectionsFigures.widthHistogram
    widthHistogram(D_cell, dD_cell, area_mat, name);
end

fprintf("    1. Sections Images Generation (%s) took %ds\n", name, round(toc))

% 2. Blood Volume Rate Figures
tic

[Q_t, dQ_t] = plotRadius(Q_mat, dQ_mat, t, index_start, index_end, name);
if contains(name, 'Artery')
    ToolBox.Signals.add('ArterialVolumeRate', Q_t, 'µL/min', t, 's', dQ_t);
else
    ToolBox.Signals.add('VenousVolumeRate', Q_t, 'µL/min', t, 's', dQ_t);
end

if params.json.CrossSectionsFigures.BloodFlowProfiles
    interpolatedBloodVelocityProfile(v_profiles_cell, dv_profiles_cell, sysIdx, diasIdx, numSections, name)
end

% Call for arterial analysis
graphCombined(M0_ff_video, imdilate(mask, strel('disk', params.json.PulseAnalysis.LocalBackgroundWidth)), ...
    Q_t, dQ_t, xy_barycenter, sprintf('%s_vr', name), ...
    'etiquettes_locs', [], ...
    'etiquettes_values', [], ...
    'ylabl', 'Volume Rate (µL/min)', ...
    'xlabl', 'Time (s)', ...
    'fig_title', 'Blood Volume Rate', ...
    'unit', 'µL/min', ...
    'skip', ~params.exportVideos, ...
    'Color', name, ...
    'Visible', false);

fprintf("    2. Blood Volume Rate Figures (%s) took %ds\n", name, round(toc))

% 3. Arterial Indicators
tic

if params.json.CrossSectionsFigures.strokeAndTotalVolume && ~isempty(systolesIndexes)
    strokeAndTotalVolume(Q_t, dQ_t, systolesIndexes, t, 1000, name);
end

if params.json.CrossSectionsFigures.ARIBVR
    ArterialResistivityIndex(t, Q_t, dQ_t, sysIdx, diasIdx, sprintf('BVR%s', name), 'crossSectionsAnalysis');
end

fprintf("    3. Arterial Indicators Images Generation (%s) took %ds\n", name, round(toc))

close all

end
