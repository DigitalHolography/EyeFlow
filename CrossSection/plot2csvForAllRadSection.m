function plot2csvForAllRadSection(t, Q_cell, Q_se_cell, Q_branch, Q_se_branch, Q_radius, Q_se_radius, initial)

[numC, numB] = size(Q_cell);
ToolBox = getGlobalToolBox;

T = table();

if not(isempty(t)) % if time is a variable
    T.time = t';
end

for cIdx = 1:numC

    for bIdx = 1:numB

        if ~isempty(Q_cell{cIdx, bIdx})
            T.(sprintf('AVGVolumeRate_R%d_B%d_%s', cIdx, bIdx, initial)) = Q_cell{cIdx, bIdx}';
            T.(sprintf('STDVolumeRate_R%d_B%d_%s', cIdx, bIdx, initial)) = Q_se_cell{cIdx, bIdx}';
        end

    end

end

for bIdx = 1:numB
    T.(sprintf('AVGVolumeRate_B%d_Total_%s', bIdx, initial)) = Q_branch(bIdx, :)';
    T.(sprintf('STDVolumeRate_B%d_Total_%s', bIdx, initial)) = Q_se_branch(bIdx, :)';
end

for bIdx = 1:numC
    T.(sprintf('AVGVolumeRate_R%d_Total_%s', cIdx, initial)) = Q_radius(cIdx, :)';
    T.(sprintf('STDVolumeRate_R%d_Total_%s', cIdx, initial)) = Q_se_radius(cIdx, :)';
end

writetable(T, fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'BloodVolumeRateTable', '_', initial, '.csv')));
end
