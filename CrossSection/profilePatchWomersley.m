function profilePatchWomersley(v_profiles_cell, name, locsLabel, M0_ff_img)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;

% if ~saveFigures
%     return;
% end

% Check sizes and extract numFrames from first non empty profile data in input
[rows, cols] = size(locsLabel);
assert(isequal(size(v_profiles_cell), [rows, cols]), 'Size of v_profiles_cell must match locsLabel');
ind = 0;
numFrames = 0;

while numFrames <= 0
    ind = ind + 1;
    numFrames = size(v_profiles_cell{ind}, 2);

    if ind > size(v_profiles_cell, 1) * size(v_profiles_cell, 2)
        warning("Velocity profiles cells are all empty.")
        break
    end

end

% Extract cardiac frequency and corresponding indices with a margin
cardiac_frequency = ToolBox.Cache.HeartBeatFFT; % in Hz

f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, numFrames);

[~, cardiac_idx] = min(abs(f - cardiac_frequency));

margin_ = round(0.1 * (cardiac_idx - numFrames / 2)); % +- 10 % of Heartrate
cardiac_idxs = cardiac_idx + (-margin_:margin_);
cardiac_idxs(cardiac_idxs > numFrames) = [];
cardiac_idxs(cardiac_idxs < 1) = [];

% Start plotting
tmp = v_profiles_cell{ind};
sizeProfiles = size(tmp{ind}, 2) * 2/3;

fi = figure("Visible", "off", 'Color', 'w');
imshow(M0_ff_img, []);
axis image
axis off
fi.Position = [200 200 600 600];

hold on;
title(['Womersley Profiles Overlay - ' name]);

% Parameters for Profiles size
% profWidth = 40;
profHeight = 30;

% AVG Plot
% lines_cell = cell(rows, cols);

idx = 1;

MAX_HARMONIC = params.json.exportCrossSectionResults.WomersleyMaxHarmonic;

for circleIdx = 1:rows

    for branchIdx = 1:cols

        % 1. Plot cardiac profiles

        if isempty(locsLabel{circleIdx, branchIdx}) || isempty(v_profiles_cell{circleIdx, branchIdx})
            continue;
        end

        % Get prof data
        profData = v_profiles_cell{circleIdx, branchIdx};

        if ~isequal(size(profData, 2), numFrames)
            warning('Expected v_profiles_cell{%d,%d} to be profile size, numFrames', circleIdx, branchIdx);
            continue;
        end

        % Calculate FFT of the time dependent profile
        profile_time = zeros(length(profData{1}), numFrames);

        for ff = 1:numFrames
            profile_time(:, ff) = profData{ff};
        end

        profile_ft = fftshift(fft(profile_time, [], 2), 2);

        % Calculate the complex Wom profile to plot

        profile_Wom = mean(profile_ft(:, cardiac_idxs), 2);

        profile_Wom = profile_Wom / mean(profile_Wom);

        % Compute axes center location
        pos = locsLabel{circleIdx, branchIdx}; % pos = [x, y]

        if isempty(pos) || numel(pos) ~= 2
            continue;
        end

        x = pos(1);
        y = pos(2);
        profile_Wom = profile_Wom / 5; % Normalize by 50
        n = numel(profile_Wom);
        x_axis = linspace(-sizeProfiles / 2, sizeProfiles / 2, n);

        % Plot profile
        x_plot = x + x_axis;
        y_data_r = y - real(profile_Wom) * profHeight; % Measured data (true profile)
        y_data_i = y - imag(profile_Wom) * profHeight; % Measured data (true profile)

        % Plot directly on image (no text, no axes)
        hold on;
        plot(x_plot, y_data_r, 'b', 'LineWidth', 1); % blue for real
        plot(x_plot, y_data_i, 'r', 'LineWidth', 1); % red for imag

        hold off;

        % 2. Fit cardiac profiles

        % TODO: temp fix for a single harmonic
        % womersley_results(circleIdx, branchIdx, :) = WomersleyNumberEstimation(profile_time, cardiac_frequency, name, idx, circleIdx, branchIdx);

        % Somehow safer than previous
        % temp_results = WomersleyNumberEstimation(profile_time, cardiac_frequency, name, idx, circleIdx, branchIdx);
        % reshaped_results = reshape(temp_results, 1, 1, []);
        % womersley_results(circleIdx, branchIdx, 1:numel(reshaped_results)) = reshaped_results;
        womersley_results(circleIdx, branchIdx) = WomersleyNumberEstimation(profile_time, cardiac_frequency, name, idx, circleIdx, branchIdx);

        % addStructToExtra(ToolBox.Output.Extra, "Womersley", womersley_results(idx)(1));

        ToolBox.Cache.WomersleyOut{circleIdx,branchIdx} = womersley_results(circleIdx, branchIdx);
        idx = idx + 1;
    end

end


saveWomersleyResults("Womersley/" + capitalize(name), womersley_results, get_unit(womersley_results));

% womersleyResultsAnalysis(womersley_results);


% Save figure
if saveFigures
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_velocities_womersley_profiles_overlay_%s.png", ToolBox.folder_name, name)));
end
close(fi);

close all;
end


function res = capitalize(str)
    res = char(str);
    res(1) = upper(res(1));
    res = string(res);
end


function saveWomersleyResults(BasePath, womersley_results, units_struct)
    arguments
        BasePath string
        womersley_results
        units_struct
    end

    
    saveWomersleyResults_handle(BasePath + "/M", expandStructField(womersley_results, "segments_metrics.MovingWallFixedNu"), units_struct.segments_metrics);
    saveWomersleyResults_handle(BasePath + "/R", expandStructField(womersley_results, "segments_metrics.RigidWallFixedNu"), units_struct.segments_metrics);
    saveWomersleyResults_handle(BasePath + "/M", expandStructField(womersley_results, "harmonic_metrics.MovingWallFixedNu"), units_struct.harmonic_metrics);
    saveWomersleyResults_handle(BasePath + "/R", expandStructField(womersley_results, "harmonic_metrics.RigidWallFixedNu"), units_struct.harmonic_metrics);
end

function saveWomersleyResults_handle(BasePath, womersley_cells, units_struct)
    arguments
        BasePath string
        womersley_cells cell % Input is now explicitly a Cell Array
        units_struct
    end

    [ok, msg] = checkStructConsistency(womersley_cells);
    if ~ok
        warning("saveWomersleyResults_handle: cells are not consistant! (Will skip save!)\n%s", msg);
        return;
    end

    ToolBox = getGlobalToolBox;

    if ~endsWith(BasePath, "/") && ~endsWith(BasePath, "_")
        BasePath = BasePath + "/";
    end

    % Find a valid index to get field names and types
    valid_idx = find(~cellfun(@isempty, womersley_cells), 1);
    
    if isempty(valid_idx)
        return; % The whole array is empty
    end
    
    sample_struct = womersley_cells{valid_idx};
    if ~isstruct(sample_struct)
        warning('Expected cell array of structs, found %s. Skipping %s.', class(sample_struct), BasePath);
        return;
    end

    field_names = fieldnames(sample_struct);
    
    for i = 1:numel(field_names)
        field       = field_names{i};
        u_struc_val = [];
        unit_field  = [];
        desc_field  = [];
        
        if ~isempty(units_struct) && isfield(units_struct, field)
            u_struc_val = units_struct.(field);
            if ~isstruct(u_struc_val) && numel(u_struc_val) >= 2
                unit_field = u_struc_val(2);
                desc_field = u_struc_val(1);
            end
        end

        % We map (Cell -> Struct -> Field) to (Cell -> FieldValue)
        % Result is an A x B x N cell array containing the values of this field.
        field_cells = extractFieldFromCells(womersley_cells, field);
        
        sample_val = field_cells{valid_idx};

        % Structs
        if isstruct(sample_val)
            saveWomersleyResults_handle(BasePath + field, field_cells, u_struc_val);
            continue;
        end

        
        empty_mask = cellfun(@isempty, field_cells);
        if any(empty_mask(:))
            if isnumeric(sample_val) && ~isreal(sample_val)
                filler = complex(NaN, NaN);
            else
                filler = NaN;
            end
            field_cells(empty_mask) = {filler};
        end

        try
            field_list = cell2mat(field_cells);
            
            fields_desc = ["circleIdx", "branchIdx"];
            
            if size(field_list, 3) > 1
                fields_desc = [fields_desc, "harmonic"];
            end

            ToolBox.Output.DimOut.add(BasePath + field, field_list, fields_desc, unit_field);
            ToolBox.Output.DimOut.add_attributes(BasePath + field, "Description", desc_field);
            
        catch ME
            warning("Skipping field '%s': Data dimensions inconsistent. (%s)", field, ME.message);
        end
    end
end

function out_cells = extractFieldFromCells(in_cells, field_name)
    % Extracts in_cells{i}.(field_name) into out_cells{i}
    % Keeps structure size identical.
    
    out_cells = cell(size(in_cells));
    mask = ~cellfun(@isempty, in_cells);
    
    out_cells(mask) = cellfun(@(s) s.(field_name), in_cells(mask), 'UniformOutput', false);
end

function idx = findFirstNonEmptyIdx(array)
    arguments
        array 
    end

    is_gap = @(s) all(cellfun(@isempty, struct2cell(s)));
    idx = find(~arrayfun(is_gap, array), 1);
end

function units_struct = get_unit(womersley_results)
    segments_metrics = struct(...
        "alpha_1",                  ["Womersley number",                            ""],                ...
        "alpha_n",                  ["Womersley number on harmonic",                ""],                ...
        "harmonic",                 ["Harmonic number",                             ""],                ...
        "Kappa_n",                  ["Condition fit",                               ""],                ...
        "residual_mag_RMS",         ["Residual magnitude RMS",                      "-"],               ...
        "residual_phase_RMS",       ["Residual phase RMS",                          "-"],               ...
        "residual_phase_RMS_msk",   ["Residual phase RMS (masked to reduce noise)", "-"],               ...
        "harmonic_SNR_dB",          ["Signal-to-noise ratio at harmonic n",         "dB"],              ...
        "fit_exitflag",             ["Reason the solver stopped (see lsqnonlin)",   "-"],               ...
        "R0",                       ["Baseline Vessel Radius",                      "m"],               ...
        "Rn",                       ["Radius harmonic (complex ?)",                 "m"],               ...
        "Cn",                       ["Drive Wall Gain",                             "m/s"],             ...
        "Dn",                       ["Moving Wall Gain",                            "m/s"],             ...
        "center",                   ["Center offset fit factor",                    ""],                ...
        "width",                    ["Scale fit factor",                            ""],                ...
        "omega_0",                  ["Fundamental angular frequency",               "rad/s"],           ...
        "omega_n",                  ["N-th harmonic angulat frequency",             "rad/s"],           ...
        "metrics",                  struct(...
            "RnR0_complex",             ["PWK ≈ D_n / C_n",                                         ""],                ...
            "Qn",                       ["Flow",                                                    "m3/s"],           ...
            "Gn",                       ["Gradient",                                                "Pa/m"],            ...
            "Kn",                       ["Complex Flow Gain",                                       ""],                ...
            "tau_n",                    ["Shear",                                                   "Pa"],              ...
            "AnA0",                     ["Area Puls.",                                              "- (m2?)"],         ...
            "nu_app",                   ["Viscosity (kinetic)",                                     "m2/s ? (Pa.s)"],   ...
            "mu_app",                   ["Viscosity (dinamic)",                                     "m2/s ? (Pa.s)"],   ...
            "H_GQ_n",                   ["Per-length impedance ratio",                              "Pa.s/m4"],         ...
            "H_GQ_n_Geonorm_abs",       ["Geo-normed impedance magnitude",                          "Pa.s/m2"],         ...
            "H_tauQ_n",                 ["Shear stress per unit pulsatile flow",                    "Pa.s/m3"],         ...
            "H_tauQ_n_Geonorm_abs",     ["Geo-normed Shear stress",                                 "Pa.s/m"],          ...
            "H_RQ_n",                   ["Dilatation per unit pulsatile flow",                      "s/m3"],            ...
            "H_RQ_n_Geonorm_abs",       ["Geo-normed Dilatation pulsatile flow",                    "s/m"],             ...
            "R_seg_n",                  ["Resistance-like contribution per unit length",            "Pa.s/m4"],         ...
            "C_seg_n",                  ["Compliance-like (inverse-compliance) per unit length",    "Pa.s/m4"],         ...
            "P_seg_diss_n",             ["Dissipative AC power per unit length",                    "W/m"],             ...
            "P_seg_store_n",            ["Reactive (stored) AC power per unit length",              "W/m"]              ...
            )                                                               ... 
    );

    harmonic_metrics = struct(...
        "RhoTau21",         ["Shear harmonic ratios |tau_2|/|tau_1|",               ""],    ...
        "RhoTau31",         ["Shear harmonic ratios |tau_3|/|tau_1|",               ""],    ...
        "DeltaPhiTau2",     ["Phase skewness of shear, PhiTau_2 - 2 * PhiTau_1",    "°"]    ...
    );

    units_struct.segments_metrics = segments_metrics;
    units_struct.harmonic_metrics = harmonic_metrics;

    idx = findFirstNonEmptyIdx(womersley_results);
    [res, diff] = diffStructs(units_struct, womersley_results(idx));
    if ~res
        warning_s("Womersley: The unit structure differ form the results!\n%s", diff);
    end
end


    
function res = expandStructField(StructArray, FieldName)
% EXPANDSTRUCTFIELD Extracts a field of size 1xN from an AxB struct array
% and returns an AxBxN cell array, preserving [] for empty entries.
%
% Usage:
%   Result = expandStructField(MyData, 'F');
    arguments
        StructArray 
        FieldName (1,1) string
    end

    [rows, cols] = size(StructArray);

    fieldPath = split(FieldName, ".");

    rawValues = num2cell(StructArray);

    for i = 1:length(fieldPath)
        currentField = fieldPath{i};
        
        % For every element in our current list, extract the sub-field.
        % If the field is missing or the parent is empty, return [].
        rawValues = cellfun(@(s) safeExtract(s, currentField), rawValues, 'UniformOutput', false);
    end

    nonEmptyIdx = find(~cellfun(@isempty, rawValues), 1);
    
    if isempty(nonEmptyIdx)
        res = cell(rows, cols, 0); 
        warning('Field path "%s" All fields are empty. Returning empty cell array.', FieldName);
        return;
    end
    
    firstItem = rawValues{nonEmptyIdx};
    N = length(firstItem);

    % Normalize data: Convert every entry into a 1xN Cell Array
    %    - If data exists: convert struct-array to cell-array
    %    - If data is empty: create a cell-array of empty brackets
    
    normalizeFn = @(x) prepareRow(x, N);
    
    processedCells = cellfun(normalizeFn, rawValues, 'UniformOutput', false);

    % vertcat stacks the (A*B) entries. Result is (A*B) x N.
    flatMatrix = vertcat(processedCells{:});
    
    % Reshape to A x B x N
    res = reshape(flatMatrix, [rows, cols, N]);

end

function val = safeExtract(s, f)
    % SAFEEXTRACT Safely gets field 'f' from struct 's'. Returns [] on failure.
    if isstruct(s) && isscalar(s) && isfield(s, f)
        val = s.(f);
    else
        val = [];
    end
end

function out = prepareRow(in, N)
    if isempty(in)
        % Fill holes with N separate []
        out = repmat({[]}, 1, N);
    else
        % Convert 1xN struct/array to 1xN cell array
        out = num2cell(in); 
    end
end


% DEBUG FUNCTIONS

    
function debug_struct_cells(in_cells, target_field)
    fprintf('Debugging field "%s"...\n', target_field);
    error_count = 0;
    
    % Use linear indexing to loop through everything
    for i = 1:numel(in_cells)
        s = in_cells{i};
        
        % Skip empty cells (your original code handles these fine)
        if isempty(s)
            continue; 
        end
        
        % Check if the structure actually has the field
        if ~isfield(s, target_field)
            error_count = error_count + 1;
            
            % Get list of fields it DOES have
            actual_fields = fieldnames(s);
            
            fprintf('------------------------------------------------\n');
            fprintf('MISSING FIELD at Index: %d\n', i);
            fprintf('  Expected: %s\n', target_field);
            
            if isempty(actual_fields)
                fprintf('  Actual: <Structure is empty struct with no fields>\n');
            else
                fprintf('  Actual fields present: %s\n', strjoin(actual_fields', ', '));
            end
            
            % Stop after 5 errors to avoid flooding the command window
            if error_count >= 5
                fprintf('... More errors exist (stopping debug output).\n');
                return;
            end
        end
    end
    
    if error_count == 0
        fprintf('Debug complete: All non-empty cells contain the field.\n');
    else
        fprintf('Debug complete: Found %d issues.\n', error_count);
    end
end

  
  