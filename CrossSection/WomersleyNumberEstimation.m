function results = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx)
    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;

    NUM_INTERP_POINTS = params.json.exportCrossSectionResults.InterpolationPoints;
    crossSectionLength = size(v_profile, 1);
    FWHM_um = 8;
    % There are several scatered like this in the files 
    PIXEL_SIZE = params.px_size; % in milimeters
    PIXEL_SIZE = PIXEL_SIZE * 1000; % in μm

    % TODO: Fix the psf kernel function
    psf_kernel = create_gaussian_psf_kernel(FWHM_um, NUM_INTERP_POINTS, crossSectionLength, PIXEL_SIZE);
    % Fit a simple PSF-convolved Parabolic/Plug model to get Geometry
    [geoParams, v_mean_interp] = fitGeometryOnMean(v_profile, psf_kernel, ToolBox);
    
    if isnan(geoParams.R0)
        warning('Geometry fit failed');
        return;
    end

    init_fit = struct(...
        'psf_kernel',   psf_kernel, ...
        'geoParams',    geoParams   ...
    );

    HARMONIC_NUMBER = params.json.exportCrossSectionResults.WomersleyMaxHarmonic;
    
    % TODO: Parfor does not seem to work with toolbox
    for i = 1:HARMONIC_NUMBER
        results(i) = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, i, init_fit, ToolBox);
    end
end


function fitParams = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, n_harmonic, init_fit, ToolBox)
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
    
    params = ToolBox.getParams;

    NUM_INTERP_POINTS = params.json.exportCrossSectionResults.InterpolationPoints;
    PIXEL_SIZE = params.px_size * 1000;
    
    % SYS_IDXS = ToolBox.Cache.sysIdx;
    % DIAS_IDXS = ToolBox.Cache.diasIdx;
    
    FFT_PADDING_FACTOR = 16;

    RHO_BLOOD = 1060; % Density of blood in kg/m^3

    PSF_KERNEL = init_fit.psf_kernel;
    
    metrics = struct(...
        "RnR0_complex",     complex(NaN, NaN),      ... % PWK ≈ D_n / C_n
        "RnR0_mag",         NaN,                    ...
        "RnR0_phase_deg",   NaN,                    ...
        "Q_n",              complex(NaN, NaN),      ... % Flow (mm3/s)
        "G_n",              complex(NaN, NaN),      ... % Gradient (Pa/m)
        "tau_n",            complex(NaN, NaN),      ... % Shear (Pa)
        "H_Zn",             complex(NaN, NaN),      ... % Impedance (Pa s/mm3/m)
        "H_Rn",             complex(NaN, NaN),      ... % Geom-Norm (Pa)
        "AnA0",             complex(NaN, NaN),      ...
        "nu_app",           NaN,                    ... % Viscosity (kinetic) (m2/s)
        "mu_app",           NaN                     ... % Viscosity (dinamic)
    );

    fitParams = struct(...
        "alpha_1",          NaN,                    ...
        "alpha_n",          NaN,                    ...
        "harmonic",         NaN,                    ...
        "K_cond",           NaN,                    ...
        "R0",               init_fit.geoParams.R0,  ...
        "Rn",               NaN,                    ...
        "Cn",               complex(NaN, NaN),      ...
        "Dn",               complex(NaN, NaN),      ...
        "center",           NaN,                    ...
        "width",            NaN,                    ...
        "omega_0",          NaN,                    ...
        "omega_n",          NaN,                    ...
        "metrics",          metrics                 ...
       );
    
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
    
    % [~, cardiac_idx] = min(abs(f - cardiac_frequency));
    % 
    % % Average over a small frequency band around the cardiac frequency for stability
    % freq_resolution = f(2) - f(1);
    % margin_hz = freq_resolution * 1.5; % Capture main lobe of the peak
    % cardiac_idxs = find(abs(f - cardiac_frequency) <= margin_hz);
    % 
    % if isempty(cardiac_idxs)
    %     warning('Cardiac frequency not found in FFT spectrum. Using closest peak.');
    %     cardiac_idxs = cardiac_idx;
    % end
    
    % ============================== [ FIT ] ============================ %
    
    % v_meas = mean(v_profile_ft(:, cardiac_idxs), 2);

    v_meas = extractHarmonicProfile(v_profile_ft, f, cardiac_frequency, n_harmonic);
    
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
    
    % psf_kernel = create_gaussian_psf_kernel(FWHM_um, NUM_INTERP_POINTS, crossSectionLength, PIXEL_SIZE);

    costFunctionHandle = @(p) costFun(p, x_coords, v_meas, n_harmonic, PSF_KERNEL);

    options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
    try
        [p_fit, ~, ~, ~, ~, ~, jacobian] = lsqnonlin(costFunctionHandle, p_init, lb, ub, options);

        % Fit Condition number
        fitParams.K_cond = sqrt(condest(jacobian'*jacobian));

        alpha_1             = p_fit(1);
        Cn_fit              = p_fit(2) + 1i * p_fit(3);
        Dn_fit              = p_fit(4) + 1i * p_fit(5);
        center_fit          = p_fit(6);
        width_fit           = p_fit(7);
    
        fitParams.alpha_1   = alpha_1;
        fitParams.alpha_n   = alpha_1 * sqrt(n_harmonic);
        fitParams.Cn        = Cn_fit;
        fitParams.Dn        = Dn_fit;
        fitParams.center    = center_fit;
        fitParams.width     = width_fit;
    
        fitParams.omega_0 = 2 * pi * cardiac_frequency;
        fitParams.omega_n = fitParams.omega_0 * n_harmonic;

        fitParams.harmonic = n_harmonic;
    
        Rn = PIXEL_SIZE * crossSectionLength / 2 * fitParams.width;
        fitParams.Rn = Rn;
    
        fitParams.metrics = calculate_symbols(fitParams, RHO_BLOOD);

    catch ME 
        warning('Womersley fit failed for %s (idx %d): %s', name, idx, ME.message);
        return;
    end

    uWom_fit = generate_moving_wall_model(p_fit, x_coords, n_harmonic, PSF_KERNEL);
    
    % [parabole_fit_systole, parabole_fit_diastole] = analyse_lumen_size(v_profile_good_idx_sav, SYS_IDXS, DIAS_IDXS);
    % 
    % estimated_width.systole = parabole_fit_systole;
    % estimated_width.diastole = parabole_fit_diastole;
    
    % ============================ [ Figures ] ========================== %
    
    hFig = figure("Visible", "off");
    hold on;

    title(sprintf('Womersley Fit for %s (idx %d) (Harmonic: %d)', name, idx, n_harmonic), 'Interpreter', 'none');
    
    p1 = plot(x_coords, real(v_meas), 'b-', 'LineWidth', 1, 'DisplayName', 'Measured Data');  % Measured Data (Real)
    plot(x_coords, imag(v_meas), 'r-', 'LineWidth', 1);    % Measured Data (Imag)

    p2 = plot(x_coords, real(uWom_fit), 'b--', 'LineWidth', 1, 'DisplayName', 'Model Fit'); % Model Fit (Real)
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

    lgd = legend([p1, p2], 'Location', 'best');
    title(lgd, 'Blue: Real, Red: Imaginary');

    
    pos = get(gca, 'Position'); 
    info_box_height = 0.15; 
    % Move the axes up and shrink its height to make room
    % [left, bottom, width, height]
    set(gca, 'Position', [pos(1), pos(2) + info_box_height, pos(3), pos(4) - info_box_height]);


    % fit_string = sprintf('α Womersley: %.2f\nCenter: %.2f\nWidth: %.2f', ...
    %                      fitParams.alpha, fitParams.center, fitParams.width);
    % annotation('textbox', [0.15 0.78 0.25 0.1], 'String', fit_string, ...
    %             'FitBoxToText', 'off', 'BackgroundColor', 'w', ...
    %             'EdgeColor', 'k', 'FontSize', 12, 'FontSize', 10);


    line1_str = sprintf('α     : %-10.2f     R_0       : %.4f', fitParams.alpha_n, init_fit.geoParams.R0);
    line2_str = sprintf('Center: %-10.2f     |R_n/R_0| : %.2f %%', fitParams.center, fitParams.metrics.RnR0_mag * 100);
    line3_str = sprintf('Width : %-10.2f     Phase(R_n): %.1f°', fitParams.width, fitParams.metrics.RnR0_phase_deg);
    
    fit_string = {line1_str, line2_str, line3_str};
    
    annotation('textbox', [0, 0, 1, info_box_height], ...
               'String', fit_string, ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 10);

    
    save_path = fullfile(ToolBox.path_png, 'Womersley');

    if ~isfolder(save_path)
        mkdir(save_path);
    end

    save_filename = fullfile(save_path, sprintf("%s_WomersleyFit_%s_idx%d_c%d_b%d_h%d.png", ToolBox.folder_name, name, idx, circleIdx, branchIdx, n_harmonic));
    
    try
        exportgraphics(hFig, save_filename, 'Resolution', 96);
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

% ========================== [ R0 CALCULATION ] ========================= %

function [geoParams, v_mean_interp] = fitGeometryOnMean(v_profile, psf_kernel, ToolBox)
    % Fit simple Poiseuille flow (1 - r^2) convolved with PSF to find Center and Width
    
    geoParams = struct(...
        'R0',           NaN, ... % R0           (Real)
        ...
        'center_norm',  NaN, ... % Shift        (Normalized)
        'width_norm',   NaN, ... % Width/R0     (Normalized)
        'DC_Amp',       NaN  ... % Amplitude    (Non-Normalized ?)
    );

    params = ToolBox.getParams;

    NUM_INTERP_POINTS = params.json.exportCrossSectionResults.InterpolationPoints;
    crossSectionLength = size(v_profile, 1);
    x_grid_normalized = linspace(-1, 1, NUM_INTERP_POINTS); 
    
    if crossSectionLength <= 1
         warning('fitGeometryFromMean: Profile too short'); 
         return; 
    end

    % TODO: Maybe not fit on interpolated values
    v_profile_interp = interp1(linspace(-1, 1, crossSectionLength), v_profile, x_grid_normalized);
    
    % PSF Kernel
    PIXEL_SIZE = params.px_size * 1000;

    v_mean = mean(v_profile_interp, 2);
    v_mean_interp = v_mean;

    % Model: Amplitude * (1 - ((x-center)/width)^2) * convol PSF

    %         [Amplitude,    Center,    Width]
    p_init  = [max(v_mean),  0,         0.8  ];

    lb      = [0,           -0.5,       0.1  ];
    ub      = [Inf,          0.5,       1.5  ];
    
    costFun = @(p) costFunDC(p, x_grid_normalized, v_mean, psf_kernel);
    options = optimoptions('lsqnonlin', 'Display', 'off');
    try
        p_fit = lsqnonlin(costFun, p_init, lb, ub, options);
    catch
        geoParams.R0 = NaN; 
        return;
    end
    
    width_normalized = p_fit(3);
    % Divided by 2 because the Normalized section is of length 2
    R0 = (width_normalized * (crossSectionLength * PIXEL_SIZE)) / 2; 
    
    geoParams.center_norm = p_fit(2);
    geoParams.width_norm = width_normalized;
    geoParams.R0 = R0;
    geoParams.DC_Amp = p_fit(1);
end

function err = costFunDC(p, x, v_meas, psf)
    amp     = p(1); 
    center  = p(2); 
    width   = p(3);

    r = (x - center) / width;
    
    profile = 1 - r .^ 2;
    profile(abs(r) > 1) = 0;
    profile = profile * amp;
    
    if ~isempty(psf)
        profile = conv(profile, psf, 'same');
    end
    
    err = profile(:) - v_meas(:);
end


% ========================== [ COST FUNCTIONS ] ========================= %

function res = costFun(p, x_coords, v_meas, n_harmonic, psf_kernel)
    res = [real(generate_moving_wall_model(p, x_coords, n_harmonic, psf_kernel) - v_meas.'); ...
           imag(generate_moving_wall_model(p, x_coords, n_harmonic, psf_kernel) - v_meas.')];
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

% ============================= [ WOMERSELY ] =========================== %

function model_profile = generate_moving_wall_model(p, x, n_harmonic, psf_kernel)
    % psf_kernel is optionnal

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

    Bn_profile(abs(r) > 1) = 0;
    Psi_n_profile(abs(r) > 1) = 0;

    if ~isempty(psf_kernel)
        Bn_profile_final = conv(Bn_profile, psf_kernel, 'same');
        Psi_n_profile_final = conv(Psi_n_profile, psf_kernel, 'same');
    else
        % If no PSF, use the ideal profiles.
        Bn_profile_final = Bn_profile;
        Psi_n_profile_final = Psi_n_profile;
    end
    
    profile = (Cn * Bn_profile_final) + (Dn * Psi_n_profile_final);
    
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

% ======================================================================= %
    
function v_meas = extractHarmonicProfile(v_profile_ft, f_vector, base_frequency, n_harmonic)
    % extractHarmonicProfile Extracts the complex velocity profile for a specific harmonic.
    %
    %   v_profile_ft   - The full, FFT-shifted Fourier spectrum of the velocity data.
    %   f_vector       - The corresponding frequency vector for the spectrum.
    %   base_frequency - The fundamental cardiac frequency (e.g., heart rate in Hz).
    %   n_harmonic     - The integer of the harmonic to extract (e.g., 1, 2, 3...).
    %
    % OUTPUT:
    %   v_meas         - The complex velocity profile (a column vector) averaged
    %                    over the frequency band of the specified harmonic.

    target_frequency = n_harmonic * base_frequency;

    freq_resolution = f_vector(2) - f_vector(1);
    margin_hz = freq_resolution * 1.5;
    
    harmonic_indices = find(abs(f_vector - target_frequency) <= margin_hz);

    if isempty(harmonic_indices)
        warning('Harmonic %d (%.2f Hz) not found in frequency band. Using closest single peak instead.', n_harmonic, target_frequency);
        [~, closest_idx] = min(abs(f_vector - target_frequency));
        harmonic_indices = closest_idx;
    end

    v_meas = mean(v_profile_ft(:, harmonic_indices), 2);
end

    
function psf_kernel = create_gaussian_psf_kernel(fwhm_um, num_points, cross_section_px, pixel_size_um)
    % Gaussian PSF kernel

    if fwhm_um <= 0
        % If FWHM is zero or negative, return an empty kernel to signify no convolution.
        psf_kernel = [];
        return;
    end

    total_width_um = cross_section_px * pixel_size_um;
    dx_um = total_width_um / num_points; % um per point

    if dx_um == 0
        error('Cannot create PSF kernel: The spatial resolution (dx_um) is zero.');
    end
    
    fwhm_in_points = fwhm_um / dx_um;

    sigma_in_points = fwhm_in_points / (2 * sqrt(2 * log(2)));

    kernel_radius = ceil(3 * sigma_in_points);
    kernel_size = 2 * kernel_radius + 1;

    % 'fspecial' creates a normalized kernel that sums to 1
    psf_kernel = fspecial('gaussian', [1, kernel_size], sigma_in_points);
end


% ============================== [ SYMBOLS ] ============================ %

function K = calculate_K_factor(lambda)
    % Calculates the flow moment factor K(alpha).
    % K(alpha) = 1 - (2*J1(lambda) / (lambda * J0(lambda)))
    
    if lambda == 0 || besselj(0, lambda) == 0
        K = complex(NaN, NaN); % Avoid division by zero
        return;
    end
    
    K = 1 - (2 * besselj(1, lambda)) / (lambda * besselj(0, lambda));
end


function metrics = calculate_symbols(fitParams, rho_blood)
    Cn          = fitParams.Cn;
    Dn          = fitParams.Dn;
    R0          = fitParams.R0 * 1e-6;
    alpha_1     = fitParams.alpha_1;
    alpha_n     = fitParams.alpha_n;
    omega_n     = fitParams.omega_n;

    Rn          = fitParams.Rn;

    lambda_n = (1i^(3/2)) * alpha_n;

    metrics = struct();
    % metrics = struct(...
    %     'RnR0_complex',     NaN, ... % PWK ≈ D_n / C_n
    %     'RnR0_mag',         NaN, ...
    %     'RnR0_phase_deg',   NaN, ...
    %     'Q_n',              NaN, ... % Flow (mm3/s)
    %     'G_n',              NaN, ... % Gradient (Pa/m)
    %     'tau_n',            NaN, ... % Shear (Pa)
    %     'H_Zn',             NaN, ... % Impedance (Pa s/mm3/m)
    %     'H_Rn',             NaN, ... % Geom-Norm (Pa)
    %     'AnA0',             NaN, ...
    %     'nu_app',           NaN, ... % Viscosity (kinetic) (m2/s)
    %     'mu_app',           NaN  ... % Viscosity (dinamic)
    % );

    % PWK (Rn/R0)
    if Cn ~= 0
        RnR0_complex = Dn / Cn;
    else
        RnR0_complex = complex(NaN, NaN);
    end

    metrics.RnR0_complex      = RnR0_complex;
    metrics.RnR0_mag          = abs(RnR0_complex);
    metrics.RnR0_phase_deg    = rad2deg(angle(RnR0_complex));

    % Flow (Qn)
    K_an = calculate_K_factor(lambda_n);
    Qn_m3_s = pi * R0 ^ 2 * Cn * K_an;
    metrics.Q_n = Qn_m3_s * 1e9; % units: mm³/s

    % Gradient (Gn)
    metrics.G_n = 1i * omega_n * rho_blood * Cn; % units: Pa/m

    % Viscosity (νapp)
    numerator = (R0 ^ 2) * fitParams.omega_0;
    denominator = alpha_1 ^ 2;
    
    if denominator > 0
        metrics.nu_app = numerator / denominator;
        metrics.mu_app = metrics.nu_app * rho_blood;
    else
        metrics.nu_app = NaN;
        metrics.mu_app = NaN;
    end

    % Shear (τn)
    if ~isnan(metrics.mu_app) && R0 > 0 && besselj(0, lambda_n) ~= 0
        metrics.tau_n = (metrics.mu_app * Cn * lambda_n * besselj(1, lambda_n)) / (R0 * besselj(0, lambda_n)); % units: Pa
    else
        metrics.tau_n = complex(NaN, NaN);
    end

    % Impedance (Ĥz,n)
    if Qn_m3_s ~= 0
        metrics.H_Zn = metrics.G_n / Qn_m3_s; % units: Pa·s/m⁴ (SI)
    else
        metrics.H_Zn = complex(NaN, NaN);
    end

    % Geom.-norm. (HR,n)
    % NOTE: The formula H_R,n = Ĥ_z,n * R₀² from the paper results in units
    % of [Pa·s·m], whereas the table specifies the units as [Pa].
    % This is a known inconsistency in the source document.
    % The calculation is implemented here as stated in the formula.
    metrics.H_Rn = metrics.H_Zn * (R0 ^ 2); % units: Pa·s·m

    % Area puls. (An/A0)
    metrics.AnA0 = 2 * metrics.RnR0_complex;

end


% +=====================================================================+ %
% |                                DEBUG                                | %
% +=====================================================================+ %

    
function visualizeDCFit(v_mean, geoParams, psf_kernel)
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
        final_model = conv(ideal_model, psf_kernel, 'same');
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

  