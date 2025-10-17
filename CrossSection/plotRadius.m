function [U_t, U_t_SE] = plotRadius(radius_U, radius_U_SE, fullTime, idx_start, idx_end, vessel_name, type_name)
% plotRadius - Plots the average flow rate per radius and its uncertainty
%   This function computes and plots the average flow rate and its
%   uncertainty over a specified time range for a given radius.

% Get global ToolBox and parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
save_figures = params.json.save_figures;
numCircles = params.json.generateCrossSectionSignals.NumberOfCircles; % Number of circles (radii)
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
dr = (r2 - r1) / numCircles;
rad = linspace(r1, r2 - dr, numCircles);

% Define color for shaded regions
Color_std = [0.7 0.7 0.7];

if strcmp(type_name, 'velocity')
    title_str = 'velocity';
    y_label = 'Velocity (mm/s)';
    unit = 'mm/s';
else
    title_str = 'flow_rate';
    y_label = 'Flow Rate (µL/min)';
    unit = 'µL/min';
end

% Compute time-averaged mean and standard deviation
U_r = squeeze(mean(radius_U(:, idx_start:idx_end), 2))'; % Mean over time
N = length(idx_start:idx_end);
U_r_SE = squeeze(sqrt(sum(radius_U_SE(:, idx_start:idx_end) .^ 2, 2)))' / N; % RMS of uncertainties

% Compute total Flow Rate over circles
U_t = squeeze(mean(radius_U, 1)); % Mean over circles
N = size(radius_U_SE, 1);
U_t_SE = squeeze(sqrt(sum(radius_U_SE .^ 2, 1))) / N; % RMS of uncertainties

% Compute statistics for the time range
mean_U = mean(U_t(idx_start:idx_end)); % Time-averaged mean
max_U = max(U_t(idx_start:idx_end)); % Maximum value in the range

if save_figures
    % Create shaded region for uncertainty
    curve1 = U_r + U_r_SE; % Upper bound
    curve2 = U_r - U_r_SE; % Lower bound
    rad2 = [rad, fliplr(rad)]; % X-values for fill
    inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

    % Plot time-averaged Flow Rate vs. radius
    figure("Visible", "off");
    fill(rad2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
    hold on;
    plot(rad, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
    plot(rad, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
    plot(rad, U_r, '-k', 'LineWidth', 2); % Mean curve
    yline(mean(U_r), '--k', 'LineWidth', 2, 'label', ...
        sprintf("%0.1f %s", mean(U_r), unit)); % Horizontal mean line

    % Format plot
    axis padded;
    axP = axis;
    axis tight;
    axT = axis;
    axis([axT(1), axT(2), 0, 1.07 * axP(4)]);

    box on
    set(gca, 'LineWidth', 2);
    set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
    ax = gca;
    ax.LineStyleOrderIndex = 1; % Reset if needed
    ax.SortMethod = 'depth'; % Try changing sorting method
    ax.Layer = 'top'; % This may help in some cases

    ylabel(y_label);
    xlabel('Radius (pixels)');

    % Export plot
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_radius_%s_mean_%s.png", ToolBox.folder_name, title_str, vessel_name)));

    % Plot radial variations of Flow Rate over time
    figure("Visible", "off");
    hold on;

    for circleIdx = 1:numCircles
        plot(fullTime, radius_U(circleIdx, :));
    end

    plot(fullTime, mean(radius_U, 1), 'k', 'LineWidth', 2);

    axis padded;
    axP = axis;
    axis tight;
    axT = axis;
    axis([axT(1), axT(2), 0, 1.07 * axP(4)]);
    box on;
    ylabel(y_label);
    xlabel('Time (s)');
    set(gca, 'PlotBoxAspectRatio', [2.5 1 1], 'LineWidth', 2);

    % Export plot
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_radius_%s_variance_%s.png", ToolBox.folder_name, title_str, vessel_name)));

    % Plot total Flow Rate over time
    figure("Visible", "off");
    curve1 = U_t + U_t_SE; % Upper bound
    curve2 = U_t - U_t_SE; % Lower bound
    ft2 = [fullTime, fliplr(fullTime)]; % X-values for fill
    inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

    hold on;
    fill(ft2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
    yline(0, 'k-', 'LineWidth', 2); % Zero line
    plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
    plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
    plot(fullTime, U_t, '-k', 'LineWidth', 2); % Mean curve
    yline(mean_U, '--k', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);

    % Mark the time range used for averaging
    plot(fullTime(idx_start), 1.07 * max_U, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
    plot(fullTime(idx_end), 1.07 * max_U, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
    plot(fullTime(idx_start:idx_end), repmat(1.07 * max_U, 1, idx_end - idx_start + 1), '-k', 'LineWidth', 2);

    % Format plot
    axis padded;
    axP = axis;
    axis tight;
    axT = axis;
    axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);

    box on
    set(gca, 'LineWidth', 2);
    set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
    ax = gca;
    ax.LineStyleOrderIndex = 1; % Reset if needed
    ax.SortMethod = 'depth'; % Try changing sorting method
    ax.Layer = 'top'; % This may help in some cases

    % Add a text label with white background at the right edge of the plot
    ax = gca;
    xPos = ax.XLim(2) * 0.8; % Right edge of the plot
    yLen = ax.YLim(2) - ax.YLim(1);
    text(xPos, mean_U + 0.1 * yLen, sprintf("%0.1f %s", mean_U, unit), ...
        'BackgroundColor', 'w', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'Margin', 1); % Small padding

    ylabel(y_label);
    xlabel('Time (s)');

    % Export plot
    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_radius_plot_%s_%s.png", ToolBox.folder_name, title_str, vessel_name)));
    exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_radius_plot_%s_%s.eps", ToolBox.folder_name, title_str, vessel_name)));
end

end
