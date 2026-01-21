function exportProfilesToH5(name, v_cell, v_safe_cell, v_profiles_cell)
% Function to save the velocities profiles 
    arguments
        name string
        v_cell cell
        v_safe_cell cell
        v_profiles_cell cell
    end

    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;
    sys_idx_list = ToolBox.Cache.sysIdxList;

    v_mat = toArray(v_cell);
    v_safe_mat = toArray(v_safe_cell);
    v_profiles_mat = toArray4D(v_profiles_cell);

    bandLimitedSignalHarmonicCount = params.json.PulseAnalysis.BandLimitedSignalHarmonicCount;

    [velocitySignalPerBeatPerSegment_whole, velocitySignalPerBeatPerSegmentFFT_whole, velocitySignalPerBeatPerSegmentBandLimited_whole] = perBeatSignalAnalysisMat(v_safe_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    velocitySignalPerBeatPerSegment_whole            = mat2cell4D_shape(velocitySignalPerBeatPerSegment_whole);
    velocitySignalPerBeatPerSegmentFFT_whole         = mat2cell4D_shape(velocitySignalPerBeatPerSegmentFFT_whole);
    velocitySignalPerBeatPerSegmentBandLimited_whole = mat2cell4D_shape(velocitySignalPerBeatPerSegmentBandLimited_whole);

    [velocitySignalPerBeatPerSegment_trunc, velocitySignalPerBeatPerSegmentFFT_trunc, velocitySignalPerBeatPerSegmentBandLimited_trunc] = perBeatSignalAnalysisMat(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    velocitySignalPerBeatPerSegment_trunc            = mat2cell4D_shape(velocitySignalPerBeatPerSegment_trunc);
    velocitySignalPerBeatPerSegmentFFT_trunc         = mat2cell4D_shape(velocitySignalPerBeatPerSegmentFFT_trunc);
    velocitySignalPerBeatPerSegmentBandLimited_trunc = mat2cell4D_shape(velocitySignalPerBeatPerSegmentBandLimited_trunc);

    % [profilePerBeatPerSegments, profilePerBeatPerSegmentsFFT, profilePerBeatPerSegmentsBandLimited] = perBeatProfileAnalysisMat(v_profiles_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    % profilePerBeatPerSegments            = mat2cell5D_shape(profilePerBeatPerSegments);
    % profilePerBeatPerSegmentsFFT         = mat2cell5D_shape(profilePerBeatPerSegmentsFFT);
    % profilePerBeatPerSegmentsBandLimited = mat2cell5D_shape(profilePerBeatPerSegmentsBandLimited);

    % ToolBox.Output.add("profilePerBeatPerSegments"            + capitalize(name), profilePerBeatPerSegments,            h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegments",            keepSize=true);
    % ToolBox.Output.add("profilePerBeatPerSegmentsFFT_abs"     + capitalize(name), abs(profilePerBeatPerSegmentsFFT),    h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegmentsFFT_abs",     keepSize=true);
    % ToolBox.Output.add("profilePerBeatPerSegmentsFFT_arg"     + capitalize(name), angle(profilePerBeatPerSegmentsFFT),  h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegmentsFFT_arg",     keepSize=true);
    % ToolBox.Output.add("profilePerBeatPerSegmentsBandLimited" + capitalize(name), profilePerBeatPerSegmentsBandLimited, h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegmentsBandLimited", keepSize=true);

    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWhole"            + capitalize(name), velocitySignalPerBeatPerSegment_whole,            h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWhole",            keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWholeFFT_abs"     + capitalize(name), abs(velocitySignalPerBeatPerSegmentFFT_whole),    h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWholeFFT_abs",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWholeFFT_arg"     + capitalize(name), angle(velocitySignalPerBeatPerSegmentFFT_whole),  h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWholeFFT_arg",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWholeBandLimited" + capitalize(name), velocitySignalPerBeatPerSegmentBandLimited_whole, h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWholeBandLimited", keepSize=true);

    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTrunc"            + capitalize(name), velocitySignalPerBeatPerSegment_trunc,            h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTrunc",            keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTruncFFT_abs"     + capitalize(name), abs(velocitySignalPerBeatPerSegmentFFT_trunc),    h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTruncFFT_abs",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTruncFFT_arg"     + capitalize(name), angle(velocitySignalPerBeatPerSegmentFFT_trunc),  h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTruncFFT_arg",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTruncBandLimited" + capitalize(name), velocitySignalPerBeatPerSegmentBandLimited_trunc, h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTruncBandLimited", keepSize=true);

    ToolBox.Output.add("velocity_trunc_seg_mean_" + name, v_mat,      h5path = capitalize(name) + "/CrossSections/velocity_trunc_seg_mean");
    ToolBox.Output.add("velocity_whole_seg_mean_" + name, v_safe_mat, h5path = capitalize(name) + "/CrossSections/velocity_whole_seg_mean");

    ToolBox.Output.add("velocity_profiles_whole_seg" + name, v_profiles_mat, h5path = capitalize(name) + "/CrossSections/velocity_profiles_whole_seg", keepSize = true);
end

function [velocitySignalPerBeat, velocitySignalPerBeatFFT, velocitySignalPerBeatBandLimited] = perBeatSignalAnalysisMat(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount)
    arguments
        v_mat,
        sys_idx_list,
        bandLimitedSignalHarmonicCount
    end
    [c_size, b_size, ~] = size(v_mat);
    velocitySignalPerBeat            = cell(c_size, b_size);
    velocitySignalPerBeatFFT         = cell(c_size, b_size);
    velocitySignalPerBeatBandLimited = cell(c_size, b_size);

    for c_idx = 1:c_size
        for b_idx = 1:b_size
            [velocitySignalPerBeat{c_idx, b_idx}, velocitySignalPerBeatFFT{c_idx, b_idx}, velocitySignalPerBeatBandLimited{c_idx, b_idx}] = perBeatSignalAnalysis(v_mat(c_idx, b_idx, :), sys_idx_list, bandLimitedSignalHarmonicCount);
        end
    end
end

function [profilePerBeat, profilePerBeatFFT, profilePerBeatBandLimited] = perBeatProfileAnalysisMat(v_profiles_4d, sys_idx_list, bandLimitedSignalHarmonicCount)
    arguments
        v_profiles_4d
        sys_idx_list
        bandLimitedSignalHarmonicCount
    end

    [c_size, b_size, ~, ~] = size(v_profiles_4d);

    profilePerBeat            = cell(c_size, b_size);
    profilePerBeatFFT         = cell(c_size, b_size);
    profilePerBeatBandLimited = cell(c_size, b_size);

    for c_idx = 1:c_size
        for b_idx = 1:b_size

            profile = squeeze(v_profiles_4d(c_idx, b_idx, :, :));

            if all(isnan(profile), 'all')
                continue
            end

            [profilePerBeat{c_idx, b_idx}, profilePerBeatFFT{c_idx, b_idx}, profilePerBeatBandLimited{c_idx, b_idx}] = perBeatProfileAnalysis(profile, sys_idx_list, bandLimitedSignalHarmonicCount);
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

function array5D = mat2cell5D_shape(input_cell)
    arguments
        input_cell cell
    end

    % grid size
    [C, B] = size(input_cell);

    % find first non-empty cell
    firstIdx = find(~cellfun(@isempty, input_cell), 1);
    if isempty(firstIdx)
        array5D = [];
        return
    end

    % infer dimensions from first valid cell
    sample = input_cell{firstIdx};   % [beats × Nfft × pixels]
    [numBeats, Nfft, numPixels] = size(sample);

    % preallocate output
    array5D = nan(C, B, numBeats, Nfft, numPixels, 'double');

    % fill array
    for c = 1:C
        for b = 1:B
            if ~isempty(input_cell{c,b})
                array5D(c, b, :, :, :) = input_cell{c,b};
            end
        end
    end
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
