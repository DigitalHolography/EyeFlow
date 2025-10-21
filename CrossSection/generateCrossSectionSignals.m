function [results] = generateCrossSectionSignals(mask, vesselName, v_RMS, M0_ff)
% generateCrossSectionSignals Perform cross-sectional analysis of retinal vessels

% Inputs:
%   - mask: Binary mask of the vessel (artery or vein)
%   - vesselName: 'artery' or 'vein'

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% Retrieve cached variables
xy_barycenter = ToolBox.Cache.xy_barycenter;
papillaDiameter = ToolBox.Cache.papillaDiameter;
sysIdx = ToolBox.Cache.sysIdx;
diasIdx = ToolBox.Cache.diasIdx;

saveFigures = params.saveFigures;
initial = upper(vesselName(1));

[numX, numY, numFrames] = size(v_RMS);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
t = ToolBox.Cache.t;
M0_ff = rescale(M0_ff);
M0_ff_img = rescale(mean(M0_ff, 3));

% 1. Mask Sectionning for all circles

% for the all circles output
tic
numCircles = params.json.generateCrossSectionSignals.NumberOfCircles;
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
dr = (r2 - r1) / numCircles;
maskSectionCircles = zeros(numX, numY, numCircles);

if strcmp(vesselName, 'artery')
    maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, sprintf('mask%s_all_sections', vesselName), mask);
else
    maskSection = createMaskSection(ToolBox, M0_ff_img, r1, r2, xy_barycenter, sprintf('mask%s_all_sections', vesselName), [], mask);
end

[labeledVessels, numBranches] = labelVesselBranches(mask, maskSection, xy_barycenter);

if saveFigures
    cmap = jet(numBranches + 1);
    imwrite(labeledVessels + 1, cmap, fullfile(ToolBox.path_png, sprintf("%s_labeledVessels_%s.png", ToolBox.folder_name, vesselName)))
end

parfor circleIdx = 1:numCircles
    r_in = r1 + (circleIdx - 1) * dr;
    r_out = r_in + dr;
    maskSectionCircles(:, :, circleIdx) = diskMask(numX, numY, r_in, r_out, center = [x_c / numX, y_c / numY]);

    % save mask image
    if saveFigures

        if strcmp(vesselName, 'artery')
            createMaskSection(ToolBox, M0_ff_img, r_in, r_out, xy_barycenter, sprintf('mask%s_section_circle_%d', vesselName, circleIdx), mask);
        else
            createMaskSection(ToolBox, M0_ff_img, r_in, r_out, xy_barycenter, sprintf('mask%s_section_circle_%d', vesselName, circleIdx), [], mask);
        end

    end

end

fprintf("    1. Mask Sectionning for all circles (%s) output took %ds\n", vesselName, round(toc))

% 2. Properties of the sections for all circles output

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

% 3. Cross-sections analysis for all circles output

tic

if ~isfolder(fullfile(ToolBox.path_png, 'vesselSegmentLumen'))
    mkdir(ToolBox.path_png, 'vesselSegmentLumen')
end

if ~isfolder(fullfile(ToolBox.path_png, 'velocityProfiles'))
    mkdir(ToolBox.path_png, 'velocityProfiles')
end

if ~isfolder(fullfile(ToolBox.path_eps, 'velocityProfiles'))
    mkdir(ToolBox.path_eps, 'velocityProfiles')
end

% Initialisation of the cells for arteries
Q_cell = cell(numCircles, numBranches); % Average flow rate
Q_SE_cell = cell(numCircles, numBranches); % Standard deviation of flow rate
v_cell = cell(numCircles, numBranches); % Top velocity
v_SE_cell = cell(numCircles, numBranches); % Standard deviation of velocity
v_profiles_cell = cell(numCircles, numBranches); % Top velocity
v_profiles_SE_cell = cell(numCircles, numBranches); % Standard deviation of velocity
A_cell = cell(numCircles, numBranches); % Cross-sectional area
D_cell = cell(numCircles, numBranches); % Cross-section width
D_SE_cell = cell(numCircles, numBranches); % Standard deviation of cross-section width
rejected_mask = zeros(numX, numY, 3); % Cross-section mask
subImg_cell = cell(numCircles, numBranches); % Sub-images of vessels
histo_v_cell = cell(numCircles, numBranches); % Histograms of vessel velocities

% Cross-Section Analysis of the arteries
parfor c_idx = 1:numCircles

    for b_idx = 1:numBranches

        if ~isempty(locsLabel{c_idx, b_idx})
            % Call crossSectionAnalysis2
            patchName = sprintf('%s%d_C%d', initial, b_idx, c_idx);
            [cross_section_results] = crossSectionAnalysis2(ToolBox, ...
                locsLabel{c_idx, b_idx}, maskLabel{c_idx, b_idx}, ...
                xy_barycenter, v_RMS, patchName, papillaDiameter);

            % Map outputs to variables
            v_cell{c_idx, b_idx} = cross_section_results.v;
            v_SE_cell{c_idx, b_idx} = cross_section_results.v_SE;
            v_profiles_cell{c_idx, b_idx} = cross_section_results.v_profiles;
            v_profiles_SE_cell{c_idx, b_idx} = cross_section_results.v_profiles_SE;
            histo_v_cell{c_idx, b_idx} = cross_section_results.v_histo;
            Q_cell{c_idx, b_idx} = cross_section_results.Q;
            Q_SE_cell{c_idx, b_idx} = cross_section_results.Q_SE;

            A_cell{c_idx, b_idx} = cross_section_results.A;
            D_cell{c_idx, b_idx} = cross_section_results.D;
            D_SE_cell{c_idx, b_idx} = cross_section_results.D_SE;

            rejected_mask = cross_section_results.rejected_masks + rejected_mask;
            subImg_cell{c_idx, b_idx} = cross_section_results.subImg_cell;
        end

    end

end

if params.json.generateCrossSectionSignals.sectionMontage && saveFigures
    sectionMontage(subImg_cell, numSections, vesselName)
end

if saveFigures
    imwrite(rejected_mask, fullfile(ToolBox.path_png, sprintf("%s_rejected_masks_%s.png", ToolBox.folder_name, vesselName)))
end

% 3.1. Transform the cell objects generated by crossSectionAnalysisAllRad
% to arrays for the following functions

[branch_Q, branch_Q_SE] = averageBranchFlux(Q_cell, Q_SE_cell);
[radius_Q, radius_Q_SE] = averageRadiusFlux(Q_cell, Q_SE_cell);

[branch_v, branch_v_SE] = averageBranchVelocity(v_cell, v_SE_cell);
[radius_v, radius_v_SE] = averageRadiusVelocity(v_cell, v_SE_cell);

% 3.2. Creates the csv files for post processing outside Eyeflow
plot2csvForAllRadSection(t, Q_cell, Q_SE_cell, branch_Q, branch_Q_SE, radius_Q, radius_Q_SE, initial)
topvel2csv(t, v_cell, v_SE_cell, initial);

% Populate the object properties
results.locsLabel = locsLabel;
results.D_cell = D_cell;
results.A_cell = A_cell;
results.maskLabel = maskLabel;
results.v_cell = v_cell;
results.v_profiles_cell = v_profiles_cell;
results.Q_cell = Q_cell;
results.branch_Q = branch_Q;
results.radius_Q = radius_Q;
results.branch_v = branch_v;
results.radius_v = radius_v;
results.histo_v_cell = histo_v_cell;
results.labeledVessels = labeledVessels;
results.rejected_mask = rejected_mask;

% Standard Errors
results.D_SE_cell = D_SE_cell;
results.v_SE_profiles_cell = v_profiles_SE_cell;
results.branch_Q_SE = branch_Q_SE;
results.radius_Q_SE = radius_Q_SE;
results.branch_v_SE = branch_v_SE;
results.radius_v_SE = radius_v_SE;

fprintf("    3. Cross-sections analysis for all circles (%s) output took %ds\n", vesselName, round(toc))

% 4. Diameter Analysis

if params.json.generateCrossSectionSignals.diameterAnalysis

    tic

    analyzeSystoleDiastole(sysIdx, diasIdx, v_RMS, locsLabel, maskLabel, ...
        numCircles, numBranches, ToolBox, initial, xy_barycenter, papillaDiameter, vesselName, numFrames);

    fprintf("    4. Diameter Analysis (%s) output took %ds\n", vesselName, round(toc))

end

end
