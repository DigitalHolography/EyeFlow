function v_profiles_cell = TimeWarpingToPeriodic(v_profiles_cell)
    arguments
        v_profiles_cell
    end

    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;

    [rows, cols] = size(v_profiles_cell);
    sysIdxList = ToolBox.Cache.sysIdxList;

    cycleLength = 2 ^ nextpow2(max(diff(ToolBox.Cache.sysIdxList)));

    for circleIdx = 1:rows
        for branchIdx = 1:cols
            if isempty(v_profiles_cell{circleIdx, branchIdx})
                continue
            end
            v_profiles_cell{circleIdx, branchIdx} = ...
                handleSeg(...
                    v_profiles_cell{circleIdx, branchIdx}, ...
                    cycleLength, ...
                    sysIdxList, ...
                    params.json.exportCrossSectionResults.TimeWarp.SignalReplication ...
                );
        end
    end
    
    % result_mat = vertcat(v_profiles_cell{1,1}{:});
    % figure; imagesc(result_mat);
end

function warped_seg = handleSeg(seg, cycleLength, sysIdxList, duplicate_factor)
    arguments
        seg
        cycleLength
        sysIdxList
        duplicate_factor
    end

    % TODO: Cut the borders, and dupplicate the signal Factor Times

    all_bounds = sysIdxList;
    cylcleCount = numel(all_bounds) - 1;

    % start_block = sysIdxList(1) - 1;
    % end_block = length(seg) - sysIdxList(end);
    % 
    % warped_seg = cell(1, (cycleLength * cylcleCount) + start_block + end_block);
    % warped_seg(1:start_block) = seg(1:start_block);
    % warped_seg(sysIdxList(end):)

    warped_seg = cell(1, cycleLength * cylcleCount);

    for i = 1:cylcleCount
        start_idx = all_bounds(i);
        end_idx   = all_bounds(i + 1) - 1;

        if end_idx < start_idx
            continue; 
        end

        curr_cycle = seg(start_idx:end_idx);
        mat = vertcat(curr_cycle{:}); 

        baseCycleLength = size(mat, 1);

        xq = linspace(1, baseCycleLength, cycleLength);
        x  = 1:baseCycleLength;
        interpolated_mat = interp1(x, mat, xq, "linear");
        
        new_cycle_cells = num2cell(interpolated_mat, 2)';
        
        out_start = ((i-1) * cycleLength) + 1;
        out_stop  = i * cycleLength;
        warped_seg(out_start:out_stop) = new_cycle_cells;
    end

    % Duplicate the signal
    warped_seg = repmat(warped_seg, 1, duplicate_factor);
end