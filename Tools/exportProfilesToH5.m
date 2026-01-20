function exportProfilesToH5(name, v_cell, v_safe_cell, v_profiles_cell, v_profiles_cropped_cell)
% Function to save the velocities profiles 
    arguments
        name string
        v_cell cell
        v_safe_cell cell
        v_profiles_cell cell
        v_profiles_cropped_cell cell
    end

    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;
    sys_idx_list = ToolBox.Cache.sysIdxList;

    v_mat = toArray(v_cell);
    v_safe_mat = toArray(v_safe_cell);

    bandLimitedSignalHarmonicCount = params.json.PulseAnalysis.BandLimitedSignalHarmonicCount;

    [velocitySignalPerBeatPerSegments, velocitySignalPerBeatPerSegmentsFFT, velocitySignalPerBeatPerSegmentsBandLimited] = perBeatAnalysisMat(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    velocitySignalPerBeatPerSegments            = mat2cell4D_shape(velocitySignalPerBeatPerSegments);
    velocitySignalPerBeatPerSegmentsFFT         = mat2cell4D_shape(velocitySignalPerBeatPerSegmentsFFT);
    velocitySignalPerBeatPerSegmentsBandLimited = mat2cell4D_shape(velocitySignalPerBeatPerSegmentsBandLimited);

    ToolBox.Output.add("velocitySignalPerBeatPerSegments"            + capitalize(name), velocitySignalPerBeatPerSegments,            h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegments",            keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentsFFT_abs"     + capitalize(name), abs(velocitySignalPerBeatPerSegmentsFFT),    h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentsFFT_abs",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentsFFT_arg"     + capitalize(name), angle(velocitySignalPerBeatPerSegmentsFFT),  h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentsFFT_arg",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentsBandLimited" + capitalize(name), velocitySignalPerBeatPerSegmentsBandLimited, h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentsBandLimited", keepSize=true);

    ToolBox.Output.add("velocity_trunc_seg_mean_" + name, v_mat,      h5path = capitalize(name) + "/CrossSections/velocity_trunc_seg_mean");
    ToolBox.Output.add("velocity_whole_seg_mean_" + name, v_safe_mat, h5path = capitalize(name) + "/CrossSections/velocity_whole_seg_mean");

    ToolBox.Output.add("velocity_profiles_whole_seg" + name, toArray4D(v_profiles_cell), h5path = capitalize(name) + "/CrossSections/velocity_profiles_whole_seg", keepSize = true);
    ToolBox.Output.add("velocity_profiles_trunc_seg" + name, toArray4D(v_profiles_cropped_cell), h5path = capitalize(name) + "/CrossSections/velocity_profiles_trunc_seg", keepSize = true);
end

function [velocitySignalPerBeat, velocitySignalPerBeatFFT, velocitySignalPerBeatBandLimited] = perBeatAnalysisMat(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount)
    arguments
        v_mat,
        sys_idx_list,
        bandLimitedSignalHarmonicCount
    end
    [c_size, b_size, ~] = size(v_mat);
    velocitySignalPerBeat = cell(c_size, b_size);
    velocitySignalPerBeatFFT = cell(c_size, b_size);
    velocitySignalPerBeatBandLimited = cell(c_size, b_size);

    for c_idx = 1:c_size
        for b_idx = 1:b_size
            [velocitySignalPerBeat{c_idx, b_idx}, velocitySignalPerBeatFFT{c_idx, b_idx}, velocitySignalPerBeatBandLimited{c_idx, b_idx}] = perBeatAnalysis(v_mat(c_idx, b_idx, :), sys_idx_list, bandLimitedSignalHarmonicCount);
        end
    end
end

function v_array = toArray(v_cell)
    [rows, cols] = size(v_cell);

    firstIdx = find(~cellfun(@isempty, v_cell), 1);
    if isempty(firstIdx)
        error("The cell array is completely empty.");
    end

    vecLen = length(v_cell{firstIdx});
    
    v_array = nan(rows, cols, vecLen, "double");
    
    for r = 1:rows
        for c = 1:cols
            if ~isempty(v_cell{r,c})
                v_array(r, c, :) = double(v_cell{r,c});
            end
        end
    end
end

function v_4d = toArray4D(v_cell)
    [rows, cols] = size(v_cell);

    firstIdx = find(~cellfun(@isempty, v_cell), 1);
    if isempty(firstIdx)
        v_4d = []; return;
    end
    
    innerCell = v_cell{firstIdx};
    numNested = numel(innerCell);
    vecLen = numel(innerCell{1});

    v_4d = nan(rows, cols, numNested, vecLen, "double");

    for i = 1:numel(v_cell)
        if ~isempty(v_cell{i})
            [r, c] = ind2sub([rows, cols], i);
            
            nestedMatrix = vertcat(v_cell{i}{:});
            v_4d(r, c, :, :) = nestedMatrix;
        end
    end
end

function array = mat2cell4D_shape(input_cell)
    arguments
        input_cell cell 
    end

    [x, y] = size(input_cell);
    [z, a] = size(input_cell{1, 1});
    array = cat(4, input_cell{:});     % z × a × (x*y)
    array = reshape(array, z, a, x, y);
    array = permute(array, [3 4 1 2]); % x × y × z × a
end

% +============================================+ %
% |                   OLD CODE                 | %
% +============================================+ %

%{
function exportProfilesToH5(name, maskLabel, v_cell, dv_cell)

% ToolBox = getGlobalToolBox;

% if nargin < 4
%     dv_cell = [];
% end

% [numCircles, numBranches] = size(maskLabel);
% % [numCircles, numBranches] = size(v_profiles_cell);

% numFrames = 0;
% i = 1;

% while numFrames <= 0
%     numFrames = size(v_cell{i}, 2);
%     i = i + 1;

%     if i > size(v_cell, 1) * size(v_cell, 2)
%         warning("Velocity profiles cells are all empty.")
%         break
%     end

% end

% L = zeros(size(maskLabel{1, 1}), 'uint8');
% idx = 1;
% % Process each circle and branch

% ProfilesTimeBranches = zeros([numFrames, numCircles, numBranches]);

% for cIdx = 1:numCircles

%     for bIdx = 1:numBranches

%         if ~isempty(v_cell{cIdx, bIdx})
%             current_v = v_cell{cIdx, bIdx};
%             cv = zeros(size(current_v{1}, 2), numFrames, "single");

%             for ff = 1:numFrames
%                 cv(:, ff) = current_v{ff};
%             end

%             if ~isempty(dv_cell)
%                 current_dv = dv_cell{cIdx, bIdx};
%                 cdv = zeros(size(current_dv{1}, 2), numFrames, "single");

%                 for ff = 1:numFrames
%                     cdv(:, ff) = current_dv{ff};
%                 end

%             else
%                 cdv = [];
%             end

%             L(maskLabel{cIdx, bIdx}) = idx;
%             ProfilesTimeBranches(:, cIdx, bIdx) = mean(cv, 1, 'omitnan');
%             ToolBox.Output.Extra.add(sprintf("Segments/%s_idx%d_c%d_b%d_VelocityProfiles", name, idx, cIdx, bIdx), cv, cdv, "mm/s");
%             idx = idx + 1;
%         end

%     end

% end

% ToolBox.Output.Extra.add(sprintf("Segments/%s_Segments_Labels", name), L);

end
%}
