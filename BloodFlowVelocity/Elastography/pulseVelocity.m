function pulseVelocity(M, ~, maskVessel, name)
ToolBox = getGlobalToolBox;
outputDir = fullfile(ToolBox.path_png, 'flexion');
params = ToolBox.getParams;

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

[L, n] = labelVesselBranches(maskVessel, ones(size(maskVessel)), ToolBox.Cache.xy_barycenter, refine = false);

if params.json.save_figures
    figure('Visible', 'off');
    imagesc(L)
    axis image; axis off;
    % Save figure
    saveas(gcf, fullfile(outputDir, ...
        sprintf("%s_%s_branches.png", ToolBox.folder_name, name)));
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
    ToolBox.Output.add("ArteryPulseWaveVelocity", PWV(idx), 'mm/s', dPWV(i));
elseif strcmp(name, 'vein')
    ToolBox.Output.add("VeinPulseWaveVelocity", PWV(idx), 'mm/s', dPWV(i));
end

close all;
end
