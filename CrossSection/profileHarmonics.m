function [DiamRatio] = profileHarmonics(v_profile, name)
% profileHarmonics creates a figure zith the different profiles esting at dif%
% INPUT:
%
%
% OUTPUT:
%   DiamRatio  - Ratio of poiseuille Diameters of average and cardiac
%   modulation.
% Create figure for static plot

% Get global ToolBox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
numFrames = size(v_profile, 2);

numInterp = params.json.exportCrossSectionResults.InterpolationPoints;
assert(size(v_profile, 1) == numInterp);

% Force the two Womersley hypothesis
%v_profile_hyp = setSymetry(v_profile);

DiamRatio = 0;

%v_profile_hyp = setBoundariesZeros(v_profile_hyp);

% Calculate Fourier transform and display
v_profile_ft = fftshift(fft(v_profile, [], 2), 2);

f = linspace(-ToolBox.fs * 1000 / ToolBox.stride / 2, ToolBox.fs * 1000 / ToolBox.stride / 2, numFrames);
% cardiac_frequency = ToolBox.Cache.HeartBeatFFT; % in Hz

harmonics = [0 ToolBox.Cache.harmonics];
harmonics_idx = zeros(1, length(harmonics));
harmonics_idx(1) = 1;

for fr = 2:length(harmonics)
    [~, idx] = min(abs(f - harmonics(fr)));
    harmonics_idx(fr) = idx;
end

% Set visualization parameters
% Color_err = [0.7 0.7 0.7];
w2w = linspace(-1, 1, numInterp);

% Create figure for static plot
figure("Visible", "off");
hold('on');

% Plot profile data
profiles = abs(v_profile_ft(:, harmonics_idx));
profiles = profiles - min(profiles, [], 1);
profiles = profiles ./ max(profiles, [], 1);
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
xline(1, '--', 'LineWidth', 2)
xline(-1, '--', 'LineWidth', 2)

baseLabels = cell(length(harmonics), 1);

for k = 1:length(harmonics)
    baseLabels{k} = sprintf("n=%d", k - 1);
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
