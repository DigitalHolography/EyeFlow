function topvel2csv(t, top_vel, top_vel_std, initial)

ToolBox = getGlobalToolBox;

%csv output of the widths
T = table();
[numC, numB] = size(top_vel);

if not(isempty(t)) % if time is a variable
    T.time = t';
end

for cIdx = 1:numC

    for bIdx = 1:numB

        if ~isempty(top_vel{cIdx, bIdx})
            T.(sprintf('Max_Vel_R%d_S%d_%s', cIdx, bIdx, initial)) = top_vel{cIdx, bIdx}';
            T.(sprintf('STD_Max_Vel_R%d_S%d_%s', cIdx, bIdx, initial)) = top_vel_std{cIdx, bIdx}';
        end

    end

end

writetable(T, fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'MaxVelocityTable', '_', initial, '.csv')));

end
