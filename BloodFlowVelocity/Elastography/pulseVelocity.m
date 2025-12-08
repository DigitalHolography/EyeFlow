function pulseVelocity(M, ~, maskVessel, name)
ToolBox = getGlobalToolBox;
outputDir = fullfile(ToolBox.path_png, 'flexion');
params = ToolBox.getParams;
saveFigures = params.saveFigures;

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

[L, n] = labelVesselBranches(maskVessel, ones(size(maskVessel)), ToolBox.Cache.xy_barycenter, refine = false);

ToolBox.Output.Extra.add(sprintf("PWV/%s_Segments_Labels", name), L);

if saveFigures
    figure('Visible', 'on');
    curcmap = jet(n+1);
    imagesc(L); colormap(curcmap);
    axis image; axis off;
    % Save figure
    saveas(gcf, fullfile(outputDir, ...
        sprintf("%s_%s_branches.png", ToolBox.folder_name, name)));
    displayBranchesWithLabels(L, save_path = fullfile(outputDir, ...
        sprintf("%s_%s_branches_names.png", ToolBox.folder_name, name)), ...
        bkgimg = ind2rgb(uint16(L), curcmap));
end

PWV = NaN(1, n);
dPWV = NaN(1, n);
scores = NaN(1, n);

for i = 1:n
    % displacementAnalysis(D, maskLongArtery);
    [PWV(i), dPWV(i), scores(i)] = pulseWaveVelocity(M, L == i, i, name);
end

[~, idx] = max(scores .* ~isnan(PWV));

if strcmp(name, 'artery')
    ToolBox.Output.add("ArteryPulseWaveVelocity", PWV(idx), 'mm/s', dPWV(idx));
elseif strcmp(name, 'vein')
    ToolBox.Output.add("VeinPulseWaveVelocity", PWV(idx), 'mm/s', dPWV(idx));
end

close all;
end
