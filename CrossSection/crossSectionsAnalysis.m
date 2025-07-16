function [Q_results] = crossSectionsAnalysis(mask, vesselName, v_RMS, M0_ff_video, xy_barycenter, papillaDiameter)

ToolBox = getGlobalToolBox;

if ~isfolder(fullfile(ToolBox.path_png, 'crossSectionsAnalysis'))
    mkdir(ToolBox.path_png, 'crossSectionsAnalysis')
    mkdir(ToolBox.path_eps, 'crossSectionsAnalysis')
end

params = ToolBox.getParams;

initial = vesselName(1);

[numX, numY, numFrames] = size(v_RMS);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
M0_ff_video = rescale(M0_ff_video);
M0_ff_img = rescale(mean(M0_ff_video, 3));

%% 1. Mask Sectionning for all circles

% for the all circles output
tic
numCircles = params.json.CrossSectionsAnalysis.NumberOfCircles;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
dr = (r2 - r1) / numCircles;
maskSectionCircles = zeros(numX, numY, numCircles);

if strcmp(vesselName, 'Artery')
    maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, sprintf('mask%s_all_sections', vesselName), mask);
else
    maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, sprintf('mask%s_all_sections', vesselName), [], mask);
end

[labeledVessels, numBranches] = labelVesselBranches(mask, maskSection, xy_barycenter);
cmap = jet(numBranches + 1);
imwrite(labeledVessels + 1, cmap, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_labeledVessels_%s.png", ToolBox.folder_name, vesselName)))

parfor circleIdx = 1:numCircles
    r_in = r1 + (circleIdx - 1) * dr;
    r_out = r_in + dr;
    maskSectionCircles(:, :, circleIdx) = diskMask(numX, numY, r_in, r_out, center = [x_c / numX, y_c / numY]);

    % save mask image
    if strcmp(vesselName, 'Artery')
        createMaskSection(ToolBox, M0_ff_img, r_in, r_out, xy_barycenter, sprintf('mask%s_section_circle_%d', vesselName, circleIdx), mask);
    else
        createMaskSection(ToolBox, M0_ff_img, r_in, r_out, xy_barycenter, sprintf('mask%s_section_circle_%d', vesselName, circleIdx), [], mask);
    end

end

fprintf("    1. Mask Sectionning for all circles (%s) output took %ds\n", vesselName, round(toc))

%% 2. Properties of the sections for all circles output

tic

locsLabel = cell(numCircles, numBranches);
maskLabel = cell(numCircles, numBranches);
numSections = zeros(1, numCircles);

parfor circleIdx = 1:numCircles
    maskSection = logical(maskSectionCircles(:, :, circleIdx) .* mask);
    s = regionprops(maskSection, 'centroid');
    numSections(circleIdx) = size(round(cat(1, s.Centroid)), 1);
end

parfor circleIdx = 1:numCircles

    for i = 1:numBranches
        maskSection = logical(maskSectionCircles(:, :, circleIdx) .* (labeledVessels == i));
        s = regionprops(maskSection, 'centroid');
        locsLabel{circleIdx, i} = round(cat(1, s.Centroid));
        maskLabel{circleIdx, i} = maskSection;
    end

end

fprintf("    2. Initialisation of the sections for all circles (%s) took %ds\n", vesselName, round(toc))

%% 3. Cross-sections analysis for all circles output

tic

if ~isfolder(fullfile(ToolBox.path_png, 'crossSectionsAnalysis', 'crossSection'))
    mkdir(fullfile(ToolBox.path_png, 'crossSectionsAnalysis'), 'crossSection')
end

if ~isfolder(fullfile(ToolBox.path_png, 'crossSectionsAnalysis', 'profiles'))
    mkdir(fullfile(ToolBox.path_png, 'crossSectionsAnalysis'), 'profiles')
    mkdir(fullfile(ToolBox.path_eps, 'crossSectionsAnalysis'), 'profiles')
end

% Initialisation of the cells for arteries
Q_cell = cell(numCircles, numBranches); % Average volume rate
Q_se_cell = cell(numCircles, numBranches); % Standard deviation of volume rate
v_cell = cell(numCircles, numBranches); % Top velocity
v_se_cell = cell(numCircles, numBranches); % Standard deviation of velocity
v_profiles_cell = cell(numCircles, numBranches); % Top velocity
v_profiles_se_cell = cell(numCircles, numBranches); % Standard deviation of velocity
A_cell = cell(numCircles, numBranches); % Cross-sectional area
D_cell = cell(numCircles, numBranches); % Cross-section width
D_se_cell = cell(numCircles, numBranches); % Standard deviation of cross-section width
rejected_mask = zeros(numX, numY, 3); % Cross-section mask
subImg_cell = cell(numCircles, numBranches); % Sub-images of vessels
histo_v_cell = cell(numCircles, numBranches); % Histograms of vessel velocities

% Cross-Section Analysis of the arteries
parfor c_idx = 1:numCircles

    for b_idx = 1:numBranches

        if ~isempty(locsLabel{c_idx, b_idx})
            % Call crossSectionAnalysis2
            patchName = sprintf('%s%d_C%d', initial, b_idx, c_idx);
            [results] = crossSectionAnalysis2(ToolBox, locsLabel{c_idx, b_idx}, maskLabel{c_idx, b_idx}, v_RMS, patchName, papillaDiameter);

            % Map outputs to variables
            v_cell{c_idx, b_idx} = results.v;
            v_se_cell{c_idx, b_idx} = results.v_se;
            v_profiles_cell{c_idx, b_idx} = results.v_profiles;
            v_profiles_se_cell{c_idx, b_idx} = results.v_profiles_se;
            histo_v_cell{c_idx, b_idx} = results.v_histo;
            Q_cell{c_idx, b_idx} = results.Q;
            Q_se_cell{c_idx, b_idx} = results.Q_se;

            A_cell{c_idx, b_idx} = results.A;
            D_cell{c_idx, b_idx} = results.D;
            D_se_cell{c_idx, b_idx} = results.D_se;

            rejected_mask = results.rejected_masks + rejected_mask;
            subImg_cell{c_idx, b_idx} = results.subImg_cell;
        end

    end

end

if params.json.CrossSectionsAnalysis.sectionMontage
    sectionMontage(subImg_cell, numSections, vesselName)
end

imwrite(rejected_mask, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_rejected_masks_%s.png", ToolBox.folder_name, vesselName)))

% 3.1. Transform the cell objects generated by crossSectionAnalysisAllRad
% to arrays for the following functions
% [area_mat, Q_mat, dQ_mat] = reshapeSections(numFrames, numSections, A_cell, Q_cell, Q_se_cell);
[branchQ, branchQSE] = averageBranchFlux(Q_cell, Q_se_cell);
[radiusQ, radiusQSE] = averageRadiusFlux(Q_cell, Q_se_cell);

% 3.2. Creates the csv files for post processing outside Eyeflow
plot2csvForAllRadSection(t, Q_cell, Q_se_cell, branchQ, branchQSE, radiusQ, radiusQSE, initial)
topvel2csv(t, v_cell, v_se_cell, initial);

Q_results.locsLabel = locsLabel;
Q_results.v_cell = v_cell;
Q_results.v_profiles_cell = v_profiles_cell;
Q_results.dv_profiles_cell = v_profiles_se_cell;
Q_results.D_cell = D_cell;
Q_results.dD_cell = D_se_cell;
Q_results.A_cell = A_cell;
Q_results.maskLabel = maskLabel;
Q_results.Q_cell = Q_cell;
Q_results.branchQ = branchQ;
Q_results.branchQSE = branchQSE;
Q_results.radiusQ = radiusQ;
Q_results.radiusQSE = radiusQSE;
Q_results.labeledVessels = labeledVessels;
Q_results.histo_v_cell = histo_v_cell;

fprintf("    3. Cross-sections analysis for all circles (%s) output took %ds\n", vesselName, round(toc))

end
