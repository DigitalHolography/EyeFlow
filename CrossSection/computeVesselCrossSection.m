function [D, dD, A, dA, c1, c2, rsquare] = computeVesselCrossSection(subImg, figName, ToolBox)

% Parameters
params = ToolBox.getParams;
px_size = params.px_size;

% Compute velocity profile
profile = mean(subImg, 1, 'omitnan');
profile(isnan(profile)) = 0;
L = length(profile);

central_range = find(profile > 0.1 * max(profile));
centt = mean(central_range);

r_range = (central_range - centt) * px_size;
[p1, p2, p3, rsquare, p1_err, p2_err, p3_err] = customPoly2Fit(r_range', profile(central_range)');
[r1, r2, r1_err, r2_err] = customPoly2Roots(p1, p2, p3, p1_err, p2_err, p3_err);

c1 = max(round(centt + (r1 / px_size)), 1);
c2 = min(round(centt + (r2 / px_size)), L);

% Determine cross-section width
D = abs(r1 - r2) / px_size;
dD = sqrt(r1_err ^ 2 + r2_err ^ 2) / px_size;

if (D > sqrt(2) * L) || (rsquare < 0.6)
    D = NaN;
    dD = NaN;
end

% Compute cross-sectional area
A = pi * (D * px_size / 2) ^ 2;
dA = pi * (px_size / 2) ^ 2 * sqrt(dD ^ 4 + 2 * dD ^ 2 * D ^ 2);

% Calculate x-axis values (position in µm)
r_ = ((1:L) - centt) * px_size * 1000;

% Calculate standard deviation and confidence interval
dprofile = std(subImg, [], 1, 'omitnan');
dprofile(isnan(dprofile)) = 0;
curve1 = profile + dprofile;
curve2 = profile - dprofile;

% Create figure
f = figAspect;

% Plot confidence interval
Color_std = [0.7, 0.7, 0.7]; % Gray color for confidence interval
fill([r_, fliplr(r_)], [curve1, fliplr(curve2)], Color_std, 'EdgeColor', 'none');
hold on;

% Plot upper and lower bounds of confidence interval
plot(r_, curve1, "Color", Color_std, 'LineWidth', 2);
plot(r_, curve2, "Color", Color_std, 'LineWidth', 2);

% Plot measured data points used for fitting
plot(r_(central_range), profile(central_range), 'xk', 'MarkerSize', 10, 'LineWidth', 1.5);

% Plot zero line
yline(0, '--k', 'LineWidth', 1);

% Plot Poiseuille fit
x_fit = linspace(min(r_), max(r_), 100) / 1000; % Interpolate for smooth fit
y_fit = p1 * x_fit .^ 2 + p2 * x_fit + p3;
plot(x_fit * 1000, y_fit, 'k', 'LineWidth', 1.5);

% Plot vessel boundaries
plot(r1 * 1000, -2, 'k|', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(r2 * 1000, -2, 'k|', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(linspace(r1 * 1000, r2 * 1000, 10), repmat(-2, 10), '-k', 'LineWidth', 1.5);

% Adjust axes
axT = axis;
axis([axT(1), axT(2), - 5, 50])
box on
set(gca, 'LineWidth', 2)
set(gca, 'PlotBoxAspectRatio', [1.618 1 1])
% Add labels and title
xlabel('Position (µm)');
ylabel('Velocity (mm/s)');
title('velocity profile and laminar flow model fit');

% Save figure

exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', 'profiles', ...
    sprintf('%s_poiseuille_profile_%s.png', ToolBox.main_foldername, figName)))

% Close figure
close(f);
end
