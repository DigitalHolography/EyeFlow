function [results] = generateCrossSectionSignals(mask, vesselName, v_RMS, M0_ff, displacementField)
% generateCrossSectionSignals Perform cross-sectional analysis of retinal vessels

% Inputs:
%   - mask: Binary mask of the vessel (artery or vein)
%   - vesselName: 'artery' or 'vein'

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% Retrieve cached variables
if ~isempty(ToolBox.Cache.xy_papilla)
    xy_barycenter = ToolBox.Cache.xy_papilla;
else
    xy_barycenter = ToolBox.Cache.xy_barycenter;
end

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
slen = params.json.generateCrossSectionSignals.SegmentsLength;
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
    r_out = r_in + slen;
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
v_safe_cell = cell(numCircles, numBranches);
v_SE_cell = cell(numCircles, numBranches); % Standard deviation of velocity
v_profiles_cell = cell(numCircles, numBranches); % Top velocity
v_profiles_cropped_cell = cell(numCircles, numBranches);
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
                xy_barycenter, v_RMS, patchName);

            % Map outputs to variables
            v_cell{c_idx, b_idx} = cross_section_results.v;
            v_safe_cell{c_idx, b_idx} = cross_section_results.v_safe;
            v_SE_cell{c_idx, b_idx} = cross_section_results.v_SE;
            v_profiles_cell{c_idx, b_idx} = cross_section_results.v_profiles;
            v_profiles_cropped_cell{c_idx, b_idx} = cross_section_results.v_profiles_cropped;
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

if params.json.Preprocess.NonRigidRegisteringFlag
    % Upscale the DispField
    dispField = displacementField.field;
    dispField = squeeze(hypot(dispField(:, :, 1, :), dispField(:, :, 2, :)));
    dispField = upscaleDispField(dispField);

    Disp_cell = cell(numCircles, numBranches);
    FFT_D_cell = cell(numCircles, numBranches);

    for c_idx = 1:numCircles

        for b_idx = 1:numBranches

            if isempty(locsLabel{c_idx, b_idx})
                continue;
            end

            patchName = sprintf('%s%d_C%d', initial, b_idx, c_idx);
            res = analyzeDispField(locsLabel{c_idx, b_idx}, maskLabel{c_idx, b_idx}, xy_barycenter, dispField, patchName);

            Disp_cell{c_idx, b_idx} = res.D;
            FFT_D_cell{c_idx, b_idx} = res.fft_D;
        end

    end

    handleDispValues(vesselName, Disp_cell, FFT_D_cell);
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
results.v_safe_cell = v_safe_cell;
results.v_profiles_cell = v_profiles_cell;
results.v_profiles_cropped_cell = v_profiles_cropped_cell;
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

% +==========================================================================+ %
% |                                  HELPER                                  | %
% +==========================================================================+ %

function v_masked_resized = upscaleDispField(dispField)

arguments
    dispField
end

ToolBox = getGlobalToolBox;
params = ToolBox.params;

% 1. Setup parameters based on your logic
[numX, numY, numFrames] = size(dispField);
kInterp = params.json.Preprocess.InterpolationFactor; % Derived from 512 -> 1023 logic

% Calculate new dimensions (matches your function's logic)
out_numX = (numX - 1) * (2 ^ kInterp - 1) + numX; % Result: 1023
out_numY = (numY - 1) * (2 ^ kInterp - 1) + numY; % Result: 1023

% 2. Pre-allocate the new array
% Using 'like' preserves data type (e.g., single vs double) to save RAM
v_masked_resized = zeros(out_numX, out_numY, numFrames, 'like', dispField);

% 3. Interpolate exactly like VideoInterpolating
% Use parfor if you have the Parallel Toolbox, otherwise change to 'for'
parfor frameIdx = 1:numFrames
    % interp2 with integer k performs the specific refinement used in your ROI generation
    v_masked_resized(:, :, frameIdx) = interp2(dispField(:, :, frameIdx), kInterp);
end

end

function handleDispValues(name, Disp_cell, ~)
[c_size, b_size] = size(Disp_cell);

% disp_peak_1_idx, disp_peak_1_val, disp_peak_2_idx, disp_peak_2_val
% res_disp = nan(c_size, b_size, 4);

% Get numFrame from first non-empty elt
idx = find(~cellfun("isempty", Disp_cell), 1);

if ~isempty(idx)
    first_data = Disp_cell{idx};
else
    warning_s("handleDispValues: All cells are empty!");
    first_data = [];
end

profile_width = size(first_data, 1);
numFrames = size(first_data, 2);

res_D_x1 = nan(c_size, b_size, numFrames);
res_D_x2 = nan(c_size, b_size, numFrames);
res_fft_x1 = nan(c_size, b_size, numFrames);
res_fft_x2 = nan(c_size, b_size, numFrames);

ToolBox = getGlobalToolBox;
fN = ToolBox.Cache.fN;

res_fft_card_freq = nan(c_size, b_size, 2);

for i = 1:c_size

    for j = 1:b_size

        current_disp = Disp_cell{i, j};

        if ~isempty(current_disp)
            % If value is complex get the abs
            disp_mean = mean(current_disp, 2);
            % res_fft_disp(i, j, :) = handlePeaks(disp_mean);
            pks = handlePeaks(disp_mean);

            x1 = pks(1);
            x2 = pks(3);

            if isnan(x1) || isnan(x2)
                warning_s("handleDispValues(%i, %i): Displacement Profile Local Max is NaN", i, j);
                continue;
            end

            % Mean with 1 idx before and after x1 & x2
            res_D_x1(i, j, :) = mean(current_disp((x1 - 1):(x1 + 1), :), 1); % current_disp(x1, :);
            res_D_x2(i, j, :) = mean(current_disp((x2 - 1):(x2 + 1), :), 1); % current_disp(x2, :);

            % Calculate mean, median, std over time
            res_x1_stats = calculate_stats_onDim(res_D_x1, 3);
            res_x2_stats = calculate_stats_onDim(res_D_x2, 3);

            % Calculate Time Spectrum of displacement Field
            res_fft_x1(i, j, :) = abs(fft(res_D_x1(i, j, :) / numFrames));
            res_fft_x2(i, j, :) = abs(fft(res_D_x2(i, j, :) / numFrames));

            cur_fft_x1 = res_fft_x1(i, j, :);
            cur_fft_x2 = res_fft_x2(i, j, :);
            res_fft_card_freq(i, j, 1) = cur_fft_x1(getCardiacIdx(res_fft_x1(i, j, :), fN));
            res_fft_card_freq(i, j, 2) = cur_fft_x2(getCardiacIdx(res_fft_x2(i, j, :), fN));
        end

        % if ~isempty(FFT_D_cell{i, j})
        %     % If value is complex get the abs
        %     FFT_seg = FFT_D_cell{i, j};
        %     res_fft_disp(i, j, :) = handlePeaks(abs(FFT_seg(:, 2)));
        % end

    end

end

res_D1_x1_stats = calculate_stats_all(res_fft_card_freq(:, :, 1));
res_D1_x2_stats = calculate_stats_all(res_fft_card_freq(:, :, 2));

% Removed Empty to reshape
emptyIndex = cellfun(@isempty, Disp_cell);
Disp_cell(emptyIndex) = {NaN(profile_width, numFrames)};


displacement_array = nan(c_size, b_size, profile_width, numFrames);

for c = 1:c_size
    for b = 1:b_size
        displacement_array(c, b, :, :) = Disp_cell{c, b};
    end
end

ToolBox.Output.add("displacementField_profile_" + name, displacement_array, h5path = capitalize(name) + "/CrossSections/DisplacementField/profile", keepSize = true);
ToolBox.Output.add("displacementField_profile_D_x1_" + name, res_D_x1, h5path = capitalize(name) + "/CrossSections/DisplacementField/profile_D_x1");
ToolBox.Output.add("displacementField_profile_D_x2_" + name, res_D_x2, h5path = capitalize(name) + "/CrossSections/DisplacementField/profile_D_x2");
ToolBox.Output.add("displacementField_profile_fft_x1_" + name, res_fft_x1, h5path = capitalize(name) + "/CrossSections/DisplacementField/profile_fft_x1");
ToolBox.Output.add("displacementField_profile_fft_x2_" + name, res_fft_x2, h5path = capitalize(name) + "/CrossSections/DisplacementField/profile_fft_x2");
ToolBox.Output.add("displacementField_D1_" + name, res_fft_card_freq, h5path = capitalize(name) + "/CrossSections/DisplacementField/D1");
ToolBox.Output.add("displacementField_D_x1_stats_" + name, res_x1_stats, h5path = capitalize(name) + "/CrossSections/DisplacementField/D_x1_stats");
ToolBox.Output.add("displacementField_D_x2_stats_" + name, res_x2_stats, h5path = capitalize(name) + "/CrossSections/DisplacementField/D_x2_stats");
ToolBox.Output.add("displacementField_D1_x1_stats_" + name, res_D1_x1_stats, h5path = capitalize(name) + "/CrossSections/DisplacementField/D1_x1_stats");
ToolBox.Output.add("displacementField_D1_x2_stats_" + name, res_D1_x2_stats, h5path = capitalize(name) + "/CrossSections/DisplacementField/D1_x2_stats");

% ToolBox.Output.DimOut.add("DispField/profile_FFT_D_local_max_" + name, res_fft_disp, []);
end

function res = handlePeaks(data)
res = [NaN, NaN, NaN, NaN];
[pks, locs] = findpeaks(data, "SortStr", "descend", "NPeaks", 2);

num_found = length(pks);

if num_found >= 1
    res(1) = locs(1);
    res(2) = pks(1);
end

if num_found >= 2
    res(3) = locs(2);
    res(4) = pks(2);
end

end

function idx = getCardiacIdx(fft_data, ~)
% GETCARDIACIDX Returns the index of the dominant heart rate peak.
%
% Usage:
%   idx = getCardiacIdx(res_fft_x1, 30);
%   heart_rate_slice = res_fft_x1(:, :, idx);

ToolBox = getGlobalToolBox;

% 1. Get size
% N = size(fft_data, 3);

% 2. Create Frequency Vector
% f = (0:N-1) * (fN / N);
f = ToolBox.Cache.f;

% 3. Define Cardiac Range (0.58 Hz (35 BPM) to 3.67 Hz (220 BPM))
%    This automatically ignores the "huge left one" (DC/Breathing)
range_idxs = find(f >= 0.58 & f <= 3.67);

% 4. Spatially average the data to find the global peak
%    (Averaging height and width makes the peak distinct from noise)
avg_spectrum = squeeze(mean(mean(fft_data, 1), 2));

% 5. Find the max peak strictly within the cardiac range
[~, local_peak] = max(avg_spectrum(range_idxs));

% 6. Convert local peak back to the global index of input array
idx = range_idxs(local_peak);

end

function res = calculate_stats_onDim(x_array, dim)

arguments
    x_array, dim
end

mea = mean(x_array, dim, "omitnan");
med = median(x_array, dim, "omitnan");
std_dev = std(x_array, 0, dim, "omitnan");

res = cat(3, mea, med, std_dev);
end

function res = calculate_stats_all(data)

arguments
    data
end

mea = mean(data, "all", "omitnan");
med = median(data, "all", "omitnan");
std_dev = std(data, 1, "all", "omitnan");

res = cat(1, mea, med, std_dev);
end
