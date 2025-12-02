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
        temp_results = WomersleyNumberEstimation(profile_time, cardiac_frequency, name, idx, circleIdx, branchIdx);
        reshaped_results = reshape(temp_results, 1, 1, []);
        womersley_results(circleIdx, branchIdx, 1:numel(reshaped_results)) = reshaped_results;

        % addStructToExtra(ToolBox.Output.Extra, "Womersley", womersley_results(idx)(1));

        ToolBox.Cache.WomersleyOut{circleIdx,branchIdx} = womersley_results(circleIdx, branchIdx);
        idx = idx + 1;
    end

end


saveWomersleyResults("Womersley/" + name, womersley_results);

womersleyResultsAnalysis(womersley_results);


% Save figure
if saveFigures
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_velocities_womersley_profiles_overlay_%s.png", ToolBox.folder_name, name)));
end
close(fi);

close all;
end


    
function saveWomersleyResults(BasePath, womersley_results)
    ToolBox = getGlobalToolBox;

    if ~endsWith(BasePath, "/") && ~endsWith(BasePath, "_")
        BasePath = BasePath + "/";
    end

    input_size = size(womersley_results);

    valid_idx = findFirstNonEmptyIdx(womersley_results);
    
    if isempty(valid_idx)
        return; % The whole array is empty
    end
    
    field_names = fieldnames(womersley_results(valid_idx));

    for i = 1:numel(field_names)
        field = field_names{i};

        raw_cells = {womersley_results.(field)};
        raw_cells = reshape(raw_cells, input_size);

        sample_val = womersley_results(valid_idx).(field);

        if isstruct(sample_val)
            empty_mask = cellfun(@isempty, raw_cells);
            
            if any(empty_mask(:))
                fnames = fieldnames(sample_val);
                struct_args = [fnames, cell(numel(fnames), 1)]'; 
                dummy_struct = struct(struct_args{:});
                
                [raw_cells{empty_mask}] = deal(dummy_struct);
            end
            
            subStructs = [raw_cells{:}];
            subStructs = reshape(subStructs, input_size);
            
            saveWomersleyResults(BasePath + field, subStructs);
            continue;
        end

        empty_mask = cellfun(@isempty, raw_cells);

        if any(empty_mask(:))
            if isnumeric(sample_val) && ~isreal(sample_val)
                filler = complex(NaN, NaN);
            else
                filler = NaN;
            end
            raw_cells(empty_mask) = {filler};
        end

        try
            field_list = cell2mat(raw_cells);
            field_list = reshape(field_list, input_size);
            
            ToolBox.Output.DimOut.add(BasePath + field, field_list, ["circleIdx", "branchIdx", "harmonic"]);
            
        catch ME
            warning("Skipping field '%s': Data dimensions inconsistent or non-scalar. (%s)", field, ME.message);
        end
    end
end


function idx = findFirstNonEmptyIdx(array)
    arguments
        array 
    end

    is_gap = @(s) all(cellfun(@isempty, struct2cell(s)));
    idx = find(~arrayfun(is_gap, array), 1);
end
  