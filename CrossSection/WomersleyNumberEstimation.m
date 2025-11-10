function results = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx)
    results = [];

    % TODO: For now a constant number of harmonics (use input parameters maybe)
    HARMONIC_NUMBER = 1;

    parfor i = 1:HARMONIC_NUMBER
        results{i} = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, i);
    end
end



function fitParams = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, n_harmonic)
    % WomersleyNumberEstimation estimates the dimensionless Womersley number (alphaWom)
    % by fitting the velocity profile to a Womersley flow model.
    %
    % 1. alpha: Womersley number
    % 2. amplitude: Complex scaling factor (magnitude and phase)
    % 3. center: Center position of the vessel
    % 4. width: Effective radius of the vessel
    %
    % INPUT:
    %   v_profile         - Array (crossSectionLength x numFrames) of velocity data.
    %   cardiac_frequency - Cardiac frequency in Hz.
    %   name, idx, ...    - Identifiers for saving the output plot.
    %
    % OUTPUT:
    %   alphaWom          - Estimated Womersley number.
    %   pseudoViscosity   - Derived dynamic viscosity in Pa·s.
    %   fitParams         - A struct containing all fitted parameters.
    
    
    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;
    NUM_INTERP_POINTS = params.json.exportCrossSectionResults.InterpolationPoints;
    PIXEL_SIZE = params.px_size;
    
    % SYS_IDXS = ToolBox.Cache.sysIdx;
    % DIAS_IDXS = ToolBox.Cache.diasIdx;
    
    FFT_PADDING_FACTOR = 16;
    
    fitParams = struct('alpha', NaN, ...
                       'pseudoViscosity', NaN, ...
                       'Cn', NaN, ...
                       'Dn', NaN, ...
                       'center', NaN, ...
                       'width', NaN, ...
                       'R1R0_complex', NaN, ...
                       'R1R0_mag', NaN, ...
                       'R1R0_phase_deg', NaN);
    % estimated_width = struct('systole', [], 'diastole', []);
    
    v_profile_avg = mean(v_profile, 2);
    valid_idxs = v_profile_avg > 0;
    v_profile = v_profile(valid_idxs, :);
    % v_profile_good_idx_sav = v_profile;
    crossSectionLength = size(v_profile, 1);
    
    if crossSectionLength > 1
        v_profile = interp1(linspace(1, crossSectionLength, crossSectionLength), v_profile, linspace(1, crossSectionLength, NUM_INTERP_POINTS));
    else
        warning('Not enough valid points in the velocity profile. Skipping fit.');
        return;
    end
    
    numFrames = size(v_profile, 2);
    N_fft = numFrames * FFT_PADDING_FACTOR;
    v_profile_ft = fftshift(fft(v_profile, N_fft, 2), 2);
    
    f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, N_fft);
    
    [~, cardiac_idx] = min(abs(f - cardiac_frequency));
    % Average over a small frequency band around the cardiac frequency for stability
    freq_resolution = f(2) - f(1);
    margin_hz = freq_resolution * 1.5; % Capture main lobe of the peak
    cardiac_idxs = find(abs(f - cardiac_frequency) <= margin_hz);
    
    if isempty(cardiac_idxs)
        warning('Cardiac frequency not found in FFT spectrum. Using closest peak.');
        cardiac_idxs = cardiac_idx;
    end
    
    % ============================== [ FIT ] ==============================
    
    v_meas = mean(v_profile_ft(:, cardiac_idxs), 2);
    
    x_coords = linspace(-1, 1, NUM_INTERP_POINTS);
    
    % p = [alpha, real(Cn), imag(Cn), real(Dn), imag(Dn), center, width]
    alpha_1_init = 4;
    Cn_init_complex = mean(v_meas(abs(v_meas)>0));
    Dn_init_complex = 0; % Start with no wall motion
    center_init = 0;
    width_init = 0.8;
    
    % p_init = [alpha_1_init, real(amp_init_complex), imag(amp_init_complex), center_init, width_init];
    p_init = [alpha_1_init, real(Cn_init_complex), imag(Cn_init_complex), real(Dn_init_complex), imag(Dn_init_complex), center_init, width_init];
    
    lb = [0.1, -Inf, -Inf, -Inf, -Inf, -0.8, 0.1]; % alpha > 0, width > 0
    ub = [20,   Inf,  Inf,  Inf,  Inf,  0.8, 1.5];
    
    options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
    try
        [p_fit, ~] = lsqnonlin(@costFun, p_init, lb, ub, options);

        alphaWom = p_fit(1);
        Cn_fit = p_fit(2) + 1i * p_fit(3);
        Dn_fit = p_fit(4) + 1i * p_fit(5);
        center_fit = p_fit(6);
        width_fit = p_fit(7);
    
        fitParams.alpha = alphaWom;
        fitParams.Cn = Cn_fit;
        fitParams.Dn = Dn_fit;
        fitParams.center = center_fit;
        fitParams.width = width_fit;

        R1R0_complex = Dn_fit / Cn_fit;
        fitParams.R1R0_complex = R1R0_complex;
        fitParams.R1R0_mag = abs(R1R0_complex);
        fitParams.R1R0_phase_deg = rad2deg(angle(R1R0_complex));
    
        omega = 2 * pi * cardiac_frequency;
    
        RHO_BLOOD = 1060; % Density of blood in kg/m^3
    
        vessel_radius_meters = PIXEL_SIZE * crossSectionLength / 2 * fitParams.width;
    
        numerator = (vessel_radius_meters^2) * omega * RHO_BLOOD;
        denominator = alphaWom ^ 2;
        
        if denominator > 0
            fitParams.pseudoViscosity = numerator / denominator;
        end
    
    catch ME
        warning('Womersley fit failed for %s (idx %d): %s', name, idx, ME.message);
        return;
    end

    uWom_fit = generate_moving_wall_model(p_fit, x_coords);
    
    % [parabole_fit_systole, parabole_fit_diastole] = analyse_lumen_size(v_profile_good_idx_sav, SYS_IDXS, DIAS_IDXS);
    % 
    % estimated_width.systole = parabole_fit_systole;
    % estimated_width.diastole = parabole_fit_diastole;
    
    % ============================ [ Figures ] ============================
    
    hFig = figure("Visible", "off");
    hold on;
    title(sprintf('Womersley Fit for %s (idx %d)', name, idx), 'Interpreter', 'none');
    plot(x_coords, real(v_meas), 'b-', 'LineWidth', 1);    % Measured Data (Real)
    plot(x_coords, imag(v_meas), 'r-', 'LineWidth', 1);    % Measured Data (Imag)
    plot(x_coords, real(uWom_fit), 'b--', 'LineWidth', 1); % Model Fit (Real)
    plot(x_coords, imag(uWom_fit), 'r--', 'LineWidth', 1); % Model Fit (Imag)
    hold off;
    
    xlim([-1 1]);
    xlabel('Normalized Cross-section', 'FontSize', 14);
    ylabel('Complex Velocity (a.u.)', 'FontSize', 14);
    legend('show', 'Location', 'best', 'FontSize', 8); 
    box on;
    grid on;
    axis tight;
    set(gca, 'LineWidth', 1.5);
    
    % fit_string = sprintf('α Womersley: %.2f\nCenter: %.2f\nWidth: %.2f', ...
    %                      fitParams.alpha, fitParams.center, fitParams.width);
    % annotation('textbox', [0.15 0.78 0.25 0.1], 'String', fit_string, ...
    %             'FitBoxToText', 'off', 'BackgroundColor', 'w', ...
    %             'EdgeColor', 'k', 'FontSize', 12, 'FontSize', 10);

    fit_string = sprintf(['α: %.2f\n' ...
                          'Center: %.2f\n' ...
                          'Width: %.2f\n' ...
                          '|R_1/R_0|: %.2f %%\n' ...
                          'Phase(R_1): %.1f°'], ...
                          fitParams.alpha, ...
                          fitParams.center, ...
                          fitParams.width, ...
                          fitParams.R1R0_mag * 100, ...
                          fitParams.R1R0_phase_deg);

    annotation('textbox', [0.15 0.75 0.3 0.15], 'String', fit_string, ...
                'FitBoxToText', 'on', 'BackgroundColor', 'w', ...
                'EdgeColor', 'k', 'FontSize', 10);
    
    save_path = fullfile(ToolBox.path_png, 'Womersley');

    if ~isfolder(save_path)
        mkdir(save_path);
    end

    save_filename = fullfile(save_path, sprintf("%s_WomersleyFit_%s_idx%d_c%d_b%d.png", ToolBox.folder_name, name, idx, circleIdx, branchIdx));
    
    try
        exportgraphics(hFig, save_filename, 'Resolution', 300);
    catch export_error
        warning('Could not save figure');
    end
    
    if ~strcmpi(get(hFig, 'Visible'), 'on')
        close(hFig);
    end
end


% +=====================================================================+ %
% |                          HELPER FUNCTIONS                           | %
% +=====================================================================+ %

% ========================== [ COST FUNCTIONS ] ===========================

function res = costFun(p, n_harmonic) 
    res = [real(generate_moving_wall_model(p, x_coords, n_harmonic) - v_meas.'); ...
           imag(generate_moving_wall_model(p, x_coords, n_harmonic) - v_meas.')];
end

function res = uWom_base(alpha, r)
    res = (1 - (besselj(0, 1i^(3/2) * alpha * r) ./ besselj(0, 1i^(3/2) * alpha)));
end


function res = uWom_psi(alpha, r)
    % Formula: Psi_n = - [lambda * J1(lambda) / J0(lambda)^2] * J0(lambda*r)
    %
    % Inputs:
    %   alpha - The Womersley number (a dimensionless real number).
    %   r     - The normalized radial coordinate (r/R0), typically a vector 
    %           from 0 to 1.
    
    % it is ((-1 + i) / sqrt(2)) * alpha
    lambda_val = 1i^(3/2) * alpha;

    numerator_scalar = -lambda_val * besselj(1, lambda_val);
    denominator_scalar = besselj(0, lambda_val) ^ 2;
    complex_scalar = numerator_scalar / denominator_scalar;
    
    radial_profile = besselj(0, lambda_val .* r);
    
    res = complex_scalar .* radial_profile;
end

% ============================= [ WOMERSELY ] =============================

function model_profile = generate_moving_wall_model(p, x, n_harmonic)
    % p = [alpha, real(Cn), imag(Cn), real(Dn), imag(Dn), center, width]
    alpha_1 = p(1);
    Cn      = p(2) + 1i * p(3);
    Dn      = p(4) + 1i * p(5);
    center  = p(6);
    width   = p(7);
    
    r = (x - center) / width;

    % See paper section 4.1, the alpha can be deduced from alpha_1
    % The functions Bn abd Psin use an alpha that is dimensionless
    % (calculated, before hand)
    alpha_n = alpha_1 * sqrt(n_harmonic);
    
    Bn_profile = uWom_base(alpha_n, r);
    Psi_n_profile = uWom_psi(alpha_n, r);
    
    profile = (Cn * Bn_profile) + (Dn * Psi_n_profile);
    
    profile(abs(r) > 1) = 0;
    
    model_profile = profile;
end

function model_profile = generate_womersley_model(p, x)
    alpha         = p(1);
    amplitude     = p(2) + 1i * p(3); % Reconstruct complex amplitude
    center        = p(4);
    width         = p(5);
    
    r = (x - center) / width;
    
    profile = uWom_base(alpha, r);
    
    profile(abs(r) > 1) = 0;
    
    model_profile = amplitude * profile;gergr
end

% =========================================================================



% +=====================================================================+ %
% |                           DEBUG FUNCTIONS                           | %
% +=====================================================================+ %

% Debug function to show parabole and its fit
function show_para_and_fit(value, fit_results)
    x_data = 1:length(value);

    p = [fit_results.p1, fit_results.p2, fit_results.p3];
    
    % Create a smooth x-axis for the fitted curve
    x_fit = linspace(fit_results.root1, fit_results.root2, 200);
    
    x_fit_original = x_fit + fit_results.center_offset;

    % Calculate the corresponding y-values for the fitted curve
    y_fit = polyval(p, x_fit);
    
    figure;
    hold on;
    
    plot(x_data, value, '.', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Full Profile Data');
    plot(fit_results.central_range, value(fit_results.central_range), 'b.', 'MarkerSize', 12, 'DisplayName', 'Data Used for Fit');
    plot(x_fit_original, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Parabolic Fit');
    
    plot([fit_results.root1, fit_results.root2] + fit_results.center_offset, [0, 0], 'kx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Roots');

    title('Comparison of Average Profile and Parabolic Fit');
    xlabel('Spatial Points (e.g., Pixel Index)');
    ylabel('Average Velocity (or other metric)');
    legend('show', 'Location', 'best'); % Add a legend
    grid on;  % Add a grid for easier reading
    hold off; % Release the plot
end

% =========================================================================
% OLD IMPLEM
% =========================================================================


% function [alphaWom, pseudoViscosity] = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx)
% % WomersleyNumberEstimation estimates the dimensionless Womersley number (alphaWom)
% % by fitting the input velocity profile (v_profile) to a Womersley flow profile.
% %
% % INPUT:
% %   v_profile - Array (crossSectionLength x numFrames) containing the measured
% % velocity profile data across time.
% %   (dv - std)
% %   cardiac_frequency - Cardiac Frequency in Hz
% %
% % OUTPUT:
% %   alphaWom  - Estimated Womersley number (dimensionless), characterizing
% %               the pulsatile flow regime.
% %   pseudoViscosity  - induced dynamic viscosity
% % Create figure for static plot
% % Get global ToolBox settings
% ToolBox = getGlobalToolBox;
% params = ToolBox.getParams;
% crossSectionLength = size(v_profile,1);
% numFrames = size(v_profile, 2);
% 
% t = ToolBox.Cache.t;
% N_padding = 16;
% 
% alphaWom = NaN;
% pseudoViscosity = NaN;
% 
% numInterp = params.json.exportCrossSectionResults.InterpolationPoints;
% 
% 
% % v_profile = interp1(linspace(1,crossSectionLength,crossSectionLength),v_profile,linspace(1,crossSectionLength,numInterp));
% 
% % Trimming the 0 on the edges
% v_profile_avg  = squeeze(mean(v_profile, 2));
% idxs = v_profile_avg~=0;
% v_profile_avg = v_profile_avg(idxs);
% v_profile = v_profile(idxs, :);
% 
% % Interpolation for more points
% sectionSize = length(v_profile_avg);
% v_profile = interp1(linspace(1, sectionSize, sectionSize), v_profile, linspace(1, sectionSize, numInterp));
% 
% % Center with a weighted mean
% v_profile_avg  = squeeze(mean(v_profile, 2));
% indxCenter = sum((1:numInterp)' .* v_profile_avg) / sum(v_profile_avg);
% shiftIndx = indxCenter - (numInterp+1)/2;
% v_profile = circshift(v_profile, round(shiftIndx), 1);
% 
% v_profile_ft = fftshift(fft(v_profile, numFrames * N_padding, 2), 2);
% 
% w2w = linspace(-1, 1, numInterp);
% 
% f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, numFrames * N_padding);
% 
% [~, cardiac_idx] = min(abs(f - cardiac_frequency));
% 
% margin_ = round(0.01 * numFrames * N_padding / 2 ); % +- 10 % fs /2
% cardiac_idxs = cardiac_idx + (-margin_:margin_);
% cardiac_idxs(cardiac_idxs > numFrames * N_padding) = [];
% cardiac_idxs(cardiac_idxs < 1) = [];
% 
% v_meas = mean(v_profile_ft(:, cardiac_idxs), 2);
% % Force the two Womersley hypothesis - bad hypothesis
% % v_profile_hyp = setSymetry(v_profile);
% 
% % v_profile_hyp = setBoundariesZeros(v_profile_hyp);  - bad hypothesis
% 
% % Calculate Fourier transform and display
% 
% % figure(1655), imagesc(t, w2w,v_profile); xlabel('Time'), ylabel('Vessel Cross Section');
% 
% % figure(1786), imagesc(f, w2w,log10(abs(v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');
% 
% % figure(1917), imagesc(angle((v_profile_ft))); xlabel('Freq'), ylabel('Vessel Cross Section ');
% 
% v_norm = v_meas / mean(v_meas);
% 
% % v_norm = v_norm - 2i * imag(v_norm);
% 
% % Create figure for static plot
% hFig = figure("Visible", "off");
% hold('on');
% 
% % Plot profile data
% plot(w2w, imag(v_norm), '-', 'Color', 'r', 'LineWidth', 2); 
% plot(w2w, real(v_norm), '-', 'Color', 'b', 'LineWidth', 2);
% 
% % Fit Womersley analytical profile
% R = 1;
% r = linspace(-1, 1, numInterp) *R; % normalized radius
% 
% uWom = @(alpha, r) (1 - (besselj(0, 1i ^ (3/2) * alpha * r) ./ besselj(0, 1i ^ (3/2) * alpha)));
% 
% % regulation_window = @(r, R) max(1 - (r / R) .^ 2, 0); % increase the weight of the central values in the fitting
% 
% uWom_tofit = @(alpha) (uWom(alpha, r) / mean(uWom(alpha, r)));
% 
% % costFun = @(alpha) norm(regulation_window(r, 1) .* uWom_tofit(alpha) - v_norm'); % least square error minimization
% costFun = @(alpha) norm(uWom_tofit(alpha) - v_norm'); % least square error minimization
% 
% alpha_init = 3; % initial guess
% alphaWom = fminsearch(costFun, alpha_init);
% % alphaWom = lsqnonlin(costFun, alpha_init, [], []);
% 
% uWom_fit = uWom_tofit(alphaWom);
% plot(w2w, imag(uWom_fit), '--', 'Color', 'r', 'LineWidth', 2);
% plot(w2w, real(uWom_fit), '--', 'Color', 'b', 'LineWidth', 2);
% 
% 
% 
% % plot(w2w, regulation_window(r, R), '--', 'Color', 'k', 'LineWidth', 2);
% 
% % Finalize static plot
% xlim([-1 1]);
% 
% xlabel('lumen cross-section (a.u.)', 'FontSize', 14);
% ylabel('Velocity (mm/s)', 'FontSize', 14);
% 
% box on
% axis tight
% axis padded;
% set(gca, 'LineWidth', 2);
% set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
% ax = gca;
% ax.LineStyleOrderIndex = 1; % Reset if needed
% ax.SortMethod = 'depth'; % Try changing sorting method
% ax.Layer = 'top'; % This may help in some cases
% 
% annotation('textbox', [0.5 0.6 0.2 0.1], ...
%             'String', sprintf('alpha Womersley : %.1f  ± %.1f', alphaWom, NaN), ...
%             'FitBoxToText', 'on', ...
%             'BackgroundColor', 'w', ...
%             'EdgeColor', 'none', ...
%             'FontSize', 12);
% 
% % Export static figure
% 
% ax = gca;
% 
% if ~isfolder(fullfile(ToolBox.path_png, 'Womersley'))
%     mkdir(fullfile(ToolBox.path_png, 'Womersley'))
% end
% 
% if isvalid(ax)
%     exportgraphics(gca, fullfile(ToolBox.path_png, "Womersley", sprintf("%s_WomersleyFit_%s_idx%d_c%d_b%d.png", ToolBox.folder_name, name, idx, circleIdx, branchIdx)), 'Resolution', 300);
% else
%     warning('Current axes are not valid. Skipping export.');
% end
% 
% % Close the figure if not needed
% if ~strcmpi(get(hFig, 'Visible'), 'on')
%     close(hFig);
% end
% 
% end
