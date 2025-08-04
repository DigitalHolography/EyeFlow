function poiseuilleProfileFigure(subImg, px_size, p1, p2, p3, r1, r2, figName, ToolBox)

% Compute velocity profile
profile = mean(subImg, 1, 'omitnan');
profile(isnan(profile)) = 0;
L = length(profile);

% Find all points above 50% threshold
central_range = find(profile > 0.2 * max(profile));
centt = mean(central_range);

% Calculate x-axis values (position in µm)
r_ = ((1:L) - centt) * px_size * 1000;

% Calculate standard deviation and confidence interval
dprofile = std(subImg, [], 1, 'omitnan');
dprofile(isnan(dprofile)) = 0;
curve1 = profile + dprofile;
curve2 = profile - dprofile;

% Create figure
f = figure("Visible", "off", "Color", 'w');

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
axis tight
axT = axis;
axis([axT(1), axT(2), - 5, 40])

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

% Add labels and title
xlabel('Position (µm)');
ylabel('Velocity (mm/s)');
fontsize(gca, 14, 'points')

% Save figure

exportgraphics(gca, fullfile(ToolBox.path_png, 'velocityProfiles', ...
    sprintf('%s_poiseuille_profile_%s.png', ToolBox.folder_name, figName)))

% Close figure
close(f);

end