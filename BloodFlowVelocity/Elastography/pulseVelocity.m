function pulseVelocity(M, ~, maskVessel, name)
% PULSEVELOCITY Estimate pulse wave velocity in vessel branches
% Inputs:
%   M           - Displacement field structure
%   maskVessel  - Binary mask of the vessel
%   name        - Name of the vessel type ('artery' or 'vein')

ToolBox = getGlobalToolBox;
outputDir = fullfile(ToolBox.path_png, 'flexion');
params = ToolBox.getParams;
saveFigures = params.saveFigures;

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

[L, n] = labelVesselBranches(maskVessel, ones(size(maskVessel)), ToolBox.Cache.xy_barycenter, refine = false);
if strcmp(name, 'artery')
    ToolBox.Output.add(sprintf("%s_PWV_Segments_Labels", name), L, h5path = "Artery/FlexuralPulseWave/BranchesLabel");
elseif strcmp(name, 'vein')
    ToolBox.Output.add(sprintf("%s_PWV_Segments_Labels", name), L, h5path = "Vein/FlexuralPulseWave/BranchesLabel");
end

if saveFigures
    figure('Visible', 'off');
    curcmap = jet(n + 1);
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
Ty = NaN(1, n);

parfor i = 1:n
    % displacementAnalysis(D, maskLongArtery);
    [PWV(i), dPWV(i), Ty(i), scores(i)] = pulseWaveVelocity(M, L == i, i, name, ToolBox);
end

if saveFigures
    displayBranchesWithLabels(L, label = PWV, unit = "mm/s", save_path = fullfile(outputDir, ...
        sprintf("%s_%s_branches_PWV.png", ToolBox.folder_name, name)), ...
        bkgimg = ind2rgb(uint16(L), curcmap));
end

if strcmp(name, 'artery')
    ToolBox.Output.add("ArteryPulseWaveVelocities", PWV, 'mm/s', dPWV, h5path = '/Artery/FlexuralPulseWave/Velocities');
elseif strcmp(name, 'vein')
    ToolBox.Output.add("VeinPulseWaveVelocities", PWV, 'mm/s', dPWV, h5path = '/Vein/FlexuralPulseWave/Velocities');
end

if strcmp(name, 'artery')
    ToolBox.Output.add("ArteryPulseWaveVelocitiesWavelengths", Ty, 'mm', h5path = '/Artery/FlexuralPulseWave/Wavelengths');
elseif strcmp(name, 'vein')
    ToolBox.Output.add("VeinPulseWaveVelocitiesWavelengths", Ty, 'mm', h5path = '/Vein/FlexuralPulseWave/Wavelengths');
end

close all;
end
