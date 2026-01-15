% +=====================================================================+ %
% |                                DEBUG                                | %
% +=====================================================================+ %

function DEBUGvisualizeDCFit(v_mean, geoParams, psf_kernel)
    % VISUALIZEDCFIT Plots the mean velocity data against the fitted model.
    % It uses the geoParams structure containing the fit results.
    %
    % INPUTS:
    %   v_mean     - The averaged velocity profile vector (the data that was fit).
    %   geoParams  - The struct output from fitGeometryOnMean.
    %   psf_kernel - The PSF kernel used during the fit.

    % --- 1. Generate the normalized coordinate grid internally ---
    num_points = length(v_mean);
    x_grid = linspace(-1, 1, num_points);

    % --- 2. Unpack parameters directly from the geoParams struct ---
    amplitude = geoParams.DC_Amp;
    center = geoParams.center_norm;
    width = geoParams.width_norm;

    % --- 3. Recreate the model components for plotting ---
    r = (x_grid - center) / width;
    ideal_model = amplitude * (1 - r.^2);
    ideal_model(abs(r) > 1) = 0;

    if ~isempty(psf_kernel)
        final_model = conv(ideal_model, psf_kernel, "same");
    else
        final_model = ideal_model;
    end

    % --- 4. Create the plot ---
    figure('Name', 'DC Fit Visualization', 'Position', [100, 100, 800, 600]);
    hold on;
    plot(x_grid, v_mean, 'b.', 'MarkerSize', 12, 'DisplayName', 'Mean Velocity Data');
    plot(x_grid, ideal_model, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ideal Parabolic Model');
    plot(x_grid, final_model, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Fitted Model (PSF-Convolved)');

    % --- 5. Add annotations ---
    xline(center, 'k:', 'DisplayName', 'Center', 'LineWidth', 2);
    xline([center - width, center + width], 'm:', 'DisplayName', 'Edges (R0)', 'HandleVisibility', 'off', 'LineWidth', 2);
    xline(center - width, 'm:', 'DisplayName', 'Edges (R0)', 'LineWidth', 2); % Re-plot one for legend
    
    grid on; box on; hold off;
    xlabel('Normalized Cross-Section Coordinate');
    ylabel('Mean Velocity');
    legend('show', 'Location', 'best');
    xlim([-1, 1]);

    fit_text = sprintf('Fitted Parameters:\n  Amplitude: %.3f\n  Center: %.3f\n  Width (R0_n): %.3f', ...
                       amplitude, center, width);
    annotation('textbox', [0.15, 0.7, 0.2, 0.2], 'String', fit_text, ...
               'BackgroundColor', 'w', 'EdgeColor', 'k', 'FitBoxToText', 'on');
end