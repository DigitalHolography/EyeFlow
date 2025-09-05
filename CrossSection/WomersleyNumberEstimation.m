function [alphaWom] = WomersleyNumberEstimation(v_profile, cardiac_frequency)
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

% Set visualization parameters
Color_err = [0.7 0.7 0.7];
w2w = linspace(-1, 1, numInterp);
% Create confidence bounds
createBounds = @(v, dv) struct( ...
    'upper', v + dv, ...
    'lower', v - dv, ...
    'x', [w2w, fliplr(w2w)], ...
    'y', [v + dv, fliplr(v - dv)]);

v_profile_ft = fftshift(fft(v_profile, [], 2), 2);

figure, imagesc(v_profile); xlabel('Time'), ylabel('Vessel Cross Section');

figure, imagesc(log10(abs(v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');

figure, imagesc(angle((v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');

f = linspace(-ToolBox.fs / 2, ToolBox.fs / 2, numFrames);

[~, cardiac_idx] = min(abs(f - cardiac_frequency));

v_meas = v_profile_ft(:, cardiac_idx);

v_norm = v_meas / mean(v_meas);

bounds = createBounds(mean(v_profile,2)', ones(1,numInterp));

% Create figure for static plot
figure("Visible", "on");
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

regulation_window = exp(-10*r.^2); % increase the weight of the central values in the fitting

costFun = @(alpha) norm(regulation_window .* uWom(alpha, r) / mean(uWom(alpha, r)) - v_norm'); % least square error minimization

alpha_init = 6; % initial guess
alphaWom = fminsearch(costFun, alpha_init);

uWom_fit = uWom(alpha_init, r);
plot(w2w, imag(uWom_fit), '--', 'Color', 'r', 'LineWidth', 2);
plot(w2w, real(uWom_fit), '--', 'Color', 'b', 'LineWidth', 2);
plot(w2w, regulation_window, '--', 'Color', 'k', 'LineWidth', 2);
warning('on', 'curvefit:fit:noStartPoint');

% Finalize static plot
xlim([-1 1]);

try
    ylim([min([bounds.lower, bounds.lower]), ...
                   1.07 * max([bounds.upper, bounds.upper])]);
catch e
    disp(e)
end
xlabel('lumen cross-section (a.u.)', 'FontSize', 14);
ylabel('Velocity (mm/s)', 'FontSize', 14);

box on
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

end
