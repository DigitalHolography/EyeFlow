function [DiamRatio] = profileHarmonics(v_profile, name)
% profileHarmonics creates a figure zith the different profiles esting at dif%
% INPUT:
%  
%
% OUTPUT:
%   DiamRatio  - Ratio of poiseuille Diameters of average and cardiac
%   modulation.
% Create figure for static plot

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
numFrames = size(v_profile, 2);

numInterp = params.json.CrossSectionsFigures.InterpolationPoints;
assert(size(v_profile, 1) == numInterp);


% Force the two Womersley hypothesis
%v_profile_hyp = setSymetry(v_profile);

DiamRatio=0;

%v_profile_hyp = setBoundariesZeros(v_profile_hyp);

% Calculate Fourier transform and display
v_profile_ft = fftshift(fft(v_profile, [], 2), 2);

f = linspace(-ToolBox.fs * 1000 / ToolBox.stride / 2, ToolBox.fs * 1000 / ToolBox.stride / 2, numFrames);
cardiac_frequency = ToolBox.Outputs.HeartBeat.value/60;

harmonics = [0 ToolBox.Cache.list.harmonics];
harmonics_idx = [];
for fr = harmonics
    [~, idx] = min(abs(f - fr));
    harmonics_idx = [harmonics_idx idx];
end

% Set visualization parameters
Color_err = [0.7 0.7 0.7];
w2w = linspace(-1, 1, numInterp);

% Create figure for static plot
figure("Visible", "on");
hold('on');

% Plot profile data
profiles = abs(v_profile_ft(:,harmonics_idx));
profiles = profiles - min(profiles,[],1); 
profiles = profiles ./ max(profiles,[],1);
plot(w2w, profiles, '-', 'LineWidth', 2);


box on
axis tight
axis padded;
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases
xline(1,'--','LineWidth', 2)
xline(-1,'--','LineWidth', 2)

baseLabels = ["n=0"];
for k=1:(length(harmonics)-1)
    baseLabels = [baseLabels;sprintf("n=%d",k)];
end
legend(baseLabels)
% Export static figure

ax = gca;

if isvalid(ax)
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_MultipleHarmonics_%s.png", ToolBox.folder_name, name)), 'Resolution', 300);
else
    warning('Current axes are not valid. Skipping export.');
end

end