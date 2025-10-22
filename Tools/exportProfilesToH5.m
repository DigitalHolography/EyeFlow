function exportProfilesToH5(name,maskLabel,v_cell,dv_cell)


ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

if nargin < 4
    dv_cell = [];
end

[numCircles, numBranches] = size(maskLabel);
% [numCircles, numBranches] = size(v_profiles_cell);

numFrames = 0;
i = 1;

while numFrames <= 0
    numFrames = size(v_cell{i}, 2);
    i = i + 1;

    if i > size(v_cell, 1) * size(v_cell, 2)
        warning("Velocity profiles cells are all empty.")
        break
    end

end

L = zeros(size(maskLabel{1,1}), 'uint8');
idx = 1;
% Process each circle and branch
for cIdx = 1:numCircles

    for bIdx = 1:numBranches
        if ~isempty(v_cell{cIdx, bIdx})
            current_v = v_cell{cIdx, bIdx};
            cv = zeros(size(current_v{1},2),numFrames,"single");
            for ff = 1:numFrames
                cv(:,ff) = current_v{ff};
            end
            if ~isempty(dv_cell)
                current_dv = dv_cell{cIdx, bIdx};
                cdv = zeros(size(current_dv{1},2),numFrames,"single");
                for ff = 1:numFrames
                    cdv(:,ff) = current_dv{ff};
                end
            else
                cdv = [];
            end
            L(maskLabel{cIdx, bIdx}) = idx;
            ToolBox.Output.Extra.add(sprintf("Segments/%s_idx%d_c%d_b%d_VelocityProfiles", name, idx, cIdx, bIdx), cv, cdv, "mm/s");
            idx = idx + 1;
        end
    end
end

ToolBox.Output.Extra.add(sprintf("Segments/%s_Segments_Labels", name), L);



end