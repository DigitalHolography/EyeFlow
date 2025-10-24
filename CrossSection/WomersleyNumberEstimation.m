function [alphaWom, pseudoViscosity] = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx)
% WomersleyNumberEstimation estimates the dimensionless Womersley number (alphaWom)
% by fitting the input velocity profile (v_profile) to a Womersley flow profile.
%
% INPUT:
%   v_profile - Array (crossSectionLength x numFrames) containing the measured
% velocity profile data across time.
%   (dv - std)
%   cardiac_frequency - Cardiac Frequency in Hz
%
% OUTPUT:
%   alphaWom  - Estimated Womersley number (dimensionless), characterizing
%               the pulsatile flow regime.
%   pseudoViscosity  - induced dynamic viscosity
% Create figure for static plot

% Get global ToolBox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
crossSectionLength = size(v_profile,1);
numFrames = size(v_profile, 2);

t = ToolBox.Cache.t;
N_padding = 16;

alphaWom = NaN;
pseudoViscosity = NaN;

numInterp = params.json.exportCrossSectionResults.InterpolationPoints;


v_profile = interp1(linspace(1,crossSectionLength,crossSectionLength),v_profile,linspace(1,crossSectionLength,numInterp));


v_profile_avg  = squeeze(mean(v_profile,2));

idxs = v_profile_avg~=0;

v_profile_avg = v_profile_avg(idxs);

v_profile = v_profile(idxs,:);


sectionSize = length(v_profile_avg);

v_profile = interp1(linspace(1,sectionSize,sectionSize),v_profile,linspace(1,sectionSize,numInterp));

v_profile_avg  = squeeze(mean(v_profile,2));

indxCenter = sum((1:numInterp)'.*v_profile_avg)/ sum(v_profile_avg);

shiftIndx = indxCenter - (numInterp+1)/2;

v_profile = circshift(v_profile, round(shiftIndx),1);

v_profile_ft = fftshift(fft(v_profile, numFrames * N_padding, 2), 2);


w2w = linspace(-1, 1, numInterp);

f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, numFrames * N_padding);

[~, cardiac_idx] = min(abs(f - cardiac_frequency));

margin_ = round(0.01 * numFrames * N_padding / 2 ); % +- 10 % fs /2
cardiac_idxs = cardiac_idx + (-margin_:margin_);
cardiac_idxs(cardiac_idxs > numFrames * N_padding) = [];
cardiac_idxs(cardiac_idxs < 1) = [];

v_meas = mean(v_profile_ft(:, cardiac_idxs),2);
% Force the two Womersley hypothesis - bad hypothesis
% v_profile_hyp = setSymetry(v_profile);

% v_profile_hyp = setBoundariesZeros(v_profile_hyp);  - bad hypothesis

% Calculate Fourier transform and display

% figure(1655), imagesc(t, w2w,v_profile); xlabel('Time'), ylabel('Vessel Cross Section');

% figure(1786), imagesc(f, w2w,log10(abs(v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');

% figure(1917), imagesc(angle((v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');

v_norm = v_meas / mean(v_meas);

% v_norm = v_norm - 2i * imag(v_norm);

% Create figure for static plot
hFig = figure("Visible", "off");
hold('on');

% Plot profile data
plot(w2w, imag(v_norm), '-', 'Color', 'r', 'LineWidth', 2); 
plot(w2w, real(v_norm), '-', 'Color', 'b', 'LineWidth', 2);

% Fit Womersley analytical profile
R = 1;
r = linspace(-1, 1, numInterp) *R; % normalized radius

uWom = @(alpha, r) (1 - (besselj(0, 1i ^ (3/2) * alpha * r) ./ besselj(0, 1i ^ (3/2) * alpha)));

% regulation_window = @(r, R) max(1 - (r / R) .^ 2, 0); % increase the weight of the central values in the fitting

uWom_tofit = @(alpha) (uWom(alpha, r) / mean(uWom(alpha, r)));

% costFun = @(alpha) norm(regulation_window(r, 1) .* uWom_tofit(alpha) - v_norm'); % least square error minimization
costFun = @(alpha) norm(uWom_tofit(alpha) - v_norm'); % least square error minimization

alpha_init = 3; % initial guess
alphaWom = fminsearch(costFun, alpha_init);

uWom_fit = uWom_tofit(alphaWom);
plot(w2w, imag(uWom_fit), '--', 'Color', 'r', 'LineWidth', 2);
plot(w2w, real(uWom_fit), '--', 'Color', 'b', 'LineWidth', 2);



% plot(w2w, regulation_window(r, R), '--', 'Color', 'k', 'LineWidth', 2);

% Finalize static plot
xlim([-1 1]);

xlabel('lumen cross-section (a.u.)', 'FontSize', 14);
ylabel('Velocity (mm/s)', 'FontSize', 14);

box on
axis tight
axis padded;
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

annotation('textbox', [0.5 0.6 0.2 0.1], ...
            'String', sprintf('alpha Womersley : %.1f  ± %.1f', alphaWom, NaN), ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'none', ...
            'FontSize', 12);

% Export static figure

ax = gca;

if ~isfolder(fullfile(ToolBox.path_png, 'Womersley'))
    mkdir(fullfile(ToolBox.path_png, 'Womersley'))
end

if isvalid(ax)
    exportgraphics(gca, fullfile(ToolBox.path_png, "Womersley", sprintf("%s_WomersleyFit_%s_idx%d_c%d_b%d.png", ToolBox.folder_name, name, idx, circleIdx, branchIdx)), 'Resolution', 300);
else
    warning('Current axes are not valid. Skipping export.');
end

% Close the figure if not needed
if ~strcmpi(get(hFig, 'Visible'), 'on')
    close(hFig);
end

end
