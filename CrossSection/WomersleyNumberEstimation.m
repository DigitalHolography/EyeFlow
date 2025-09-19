function [alphaWom] = WomersleyNumberEstimation(v_profile, cardiac_frequency, name)
% WomersleyNumberEstimation estimates the dimensionless Womersley number (alphaWom)
% by fitting the input velocity profile (v_profile) to a Womersley flow profile.
%
% INPUT:
%   v_profile - Array (numInterp x numFrames) containing the measured
% velocity profile data across time.
%   (dv - std)
%   cardiac_frequency - Cardiac Frequency in Hz
%
% OUTPUT:
%   alphaWom  - Estimated Womersley number (dimensionless), characterizing
%               the pulsatile flow regime.
% Create figure for static plot

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
numFrames = size(v_profile, 2);

numInterp = params.json.CrossSectionsFigures.InterpolationPoints;
assert(size(v_profile, 1) == numInterp);
w2w = linspace(-1, 1, numInterp);


% Force the two Womersley hypothesis
v_profile_hyp = setSymetry(v_profile);

%v_profile_hyp = setBoundariesZeros(v_profile_hyp);

% Calculate Fourier transform and display
v_profile_ft = fftshift(fft(v_profile_hyp, [], 2), 2);

% figure(1655), imagesc(v_profile); xlabel('Time'), ylabel('Vessel Cross Section');

% figure(1786), imagesc(log10(abs(v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');

% figure(1917), imagesc(angle((v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');

f = linspace(-ToolBox.fs * 1000 / ToolBox.stride / 2, ToolBox.fs * 1000 / ToolBox.stride / 2, numFrames);

[~, cardiac_idx] = min(abs(f - cardiac_frequency));

v_meas = v_profile_ft(:, cardiac_idx);

v_norm = v_meas / mean(v_meas);

% Create figure for static plot
hFig = figure("Visible", "on");
hold('on');

% Plot profile data
plot(w2w, imag(v_norm), '-', 'Color', 'r', 'LineWidth', 2);
plot(w2w, real(v_norm), '-', 'Color', 'b', 'LineWidth', 2);

% Add Womersley fits
warning('off', 'curvefit:fit:noStartPoint');

% Fit Womersley analytical profile
r = linspace(-1, 1, numInterp); % normalized radius
R = 1; % normalized vessel radius


uWom = @(alpha, r) (1 - (besselj(0, 1i^(3/2)*alpha*r/R) ./ besselj(0, 1i^(3/2)*alpha)));

regulation_window = @(r,R) max(1-(r/R).^2,0); % increase the weight of the central values in the fitting

uWom_tofit = @(alpha_multi) (uWom(alpha_multi(1), r*alpha_multi(2)) / mean(uWom(alpha_multi(1), r*alpha_multi(2))));

costFun = @(alpha_multi) norm(regulation_window(r*alpha_multi(2),1) .* uWom_tofit(alpha_multi) - v_norm'); % least square error minimization

alpha_init = [3 1 0.5]; % initial guess
alphaWom = fminsearch(costFun, alpha_init);
alpha = alphaWom(1);
r_scale = alphaWom(2);
w2w_scaled = w2w .* r_scale;
uWom_fit = uWom_tofit([3 1 0.5]);
plot(w2w_scaled, imag(uWom_fit), '--', 'Color', 'r', 'LineWidth', 2);
plot(w2w_scaled, real(uWom_fit), '--', 'Color', 'b', 'LineWidth', 2);
plot(w2w, regulation_window(r,R), '--', 'Color', 'k', 'LineWidth', 2);

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

% Export static figure

ax = gca;

if isvalid(ax)
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_WomersleyFit_%s.png", ToolBox.folder_name, name)), 'Resolution', 300);
else
    warning('Current axes are not valid. Skipping export.');
end

% Close the figure if not needed
if ~strcmpi(get(hFig, 'Visible'), 'on')
    close(hFig);
end

end
