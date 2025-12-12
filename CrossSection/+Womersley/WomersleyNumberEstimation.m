function results = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, d_profile)
    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;

    NUM_INTERP_POINTS = params.json.exportCrossSectionResults.InterpolationPoints;
    crossSectionLength = size(v_profile, 1);
    FWHM_um = 8;
    % TODO: Remove all and refactor
    % There are several scatered like this in the files (If change, need oto check others)
    PIXEL_SIZE = params.px_size; % in milimeters
    PIXEL_SIZE = PIXEL_SIZE * 1000; % in μm
    RHO_BLOOD = 1060; % Density of blood in kg/m^3

    psf_kernel = create_gaussian_psf_kernel(FWHM_um, NUM_INTERP_POINTS, crossSectionLength, PIXEL_SIZE);
    % Fit a simple PSF-convolved Parabolic/Plug model to get Geometry
    [geoParams, ~] = fitGeometryOnMean(v_profile, psf_kernel, ToolBox);
    
    if isnan(geoParams.R0)
        warning("[WOMERSLEY] Geometry fit failed");
        return;
    end

    init_fit = struct(...
        'psf_kernel',   psf_kernel, ...
        'geoParams',    geoParams   ...
    );

    HARMONIC_NUMBER = params.json.exportCrossSectionResults.WomersleyMaxHarmonic;

    % TODO: Parfor does not seem to work with toolbox
    for i = 1:HARMONIC_NUMBER
        fitParams.RigidWallFixedNu(i)  = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, i, init_fit, d_profile, ToolBox, ModelType="rigid");
        fitParams.MovingWallFixedNu(i) = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, i, init_fit, d_profile, ToolBox, ModelType="moving");
        fitParams.RigidWallFreeNu(i)   = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, i, init_fit, d_profile, ToolBox, ModelType="rigid", NuType="free");
    end

    % Quality Control
    qc.qc_RigidWallFixedNu  = womersley_qc(fitParams.RigidWallFixedNu,  "rigid");
    qc.qc_MovingWallFixedNu = womersley_qc(fitParams.MovingWallFixedNu, "moving");
    qc.qc_RigidWallFreeNu   = womersley_qc(fitParams.RigidWallFreeNu,   "rigid");

    % h_metrics.RigidWallFixedNu  = calculate_harmonic_metrics(fitParams.RigidWallFixedNu);
    % h_metrics.MovingWallFixedNu = calculate_harmonic_metrics(fitParams.MovingWallFixedNu);
    % h_metrics.RigidWallFreeNu   = calculate_harmonic_metrics(fitParams.RigidWallFreeNu);

    % filtered metrics (Should be done inside Analysis Program, not Here)
    h_metrics.RigidWallFixedNu  = calculate_harmonic_metrics(fitParams.RigidWallFixedNu(qc.qc_RigidWallFixedNu.harmonic_valid));
    h_metrics.MovingWallFixedNu = calculate_harmonic_metrics(fitParams.MovingWallFixedNu(qc.qc_MovingWallFixedNu.harmonic_valid));
    h_metrics.RigidWallFreeNu   = calculate_harmonic_metrics(fitParams.RigidWallFreeNu(qc.qc_RigidWallFreeNu.harmonic_valid));

    derived.RigidWallFixedNu  = calculate_segment_derived_metrics(fitParams.RigidWallFixedNu, qc.qc_RigidWallFixedNu, init_fit.geoParams, RHO_BLOOD);
    derived.MovingWallFixedNu = calculate_segment_derived_metrics(fitParams.MovingWallFixedNu, qc.qc_MovingWallFixedNu, init_fit.geoParams, RHO_BLOOD);
    derived.RigidWallFreeNu   = calculate_segment_derived_metrics(fitParams.RigidWallFreeNu, qc.qc_RigidWallFreeNu, init_fit.geoParams, RHO_BLOOD);

    % metrics for each segments / harmonic
    results.segments_metrics = fitParams;

    % metrics for each segments (using multiple harmonics)
    results.harmonic_metrics = h_metrics;

    % QCs are both, mostly calculated for each segments / harmonics, but some 
    % are just for segments
    results.qc = qc;

    % TODO: Maybe removed derived and put inside respective results
    % Derived are metrics that are some extra metrics, they could be included 
    % inside segment_metrics and harmonic
    % Some are segments/harmonics and some not
    results.derived = derived;
end


function fitParams = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, n_harmonic, init_fit, d_profile, ToolBox, options)
    arguments
        v_profile, cardiac_frequency, name, idx, circleIdx, branchIdx, n_harmonic, init_fit, d_profile, ToolBox
        options.ModelType   string {mustBeMember(options.ModelType, ["rigid", "moving"])} = "moving"
        options.NuType      string {mustBeMember(options.NuType,    ["fixed", "free"  ])} = "fixed"
    end
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

    fitParams = struct(...
        "alpha_1",                  NaN,                            ... % Womersley number                              (-)
        "alpha_n",                  NaN,                            ... % Womersley number on harmonic                  (-)
        "harmonic",                 NaN,                            ... % Harmonic number                               (-)
        "Kappa_n",                  NaN,                            ... % Condition fit                                 (-)
        "residual_mag_RMS",         NaN,                            ... % Residual magnitude RMS                        (-)
        "residual_phase_RMS",       NaN,                            ... % Residual phase RMS                            (-)
        "residual_phase_RMS_msk",   NaN,                            ... % Residual phase RMS (masked to reduce noise)   (-)
        "harmonic_SNR_dB",          NaN,                            ... % Signal-to-noise ratio at harmonic n           (dB)
        "fit_exitflag",             NaN,                            ... % Reason the solver stopped (see lsqnonlin)     (-)
        "R0",                       init_fit.geoParams.R0,          ... % Baseline Vessel Radius                        (m)
        "Rn",                       NaN,                            ... % Radius harmonic (complex ?)                   (m)
        "Cn",                       complex(NaN, NaN),              ... % Drive Wall Gain                               (m/s)
        "Dn",                       complex(NaN, NaN),              ... % Moving Wall Gain                              (m/s)
        "center",                   NaN,                            ... % Center offset fit factor                      (-)
        "width",                    init_fit.geoParams.width_norm,  ... % Scale fit factor                              (-)
        "omega_0",                  NaN,                            ... % Fundamental angular frequency                 (rad/s)
        "omega_n",                  NaN,                            ... % N-th harmonic angulat frequency               (rad/s)
        "metrics",                  struct(...
            "RnR0_complex",             complex(NaN, NaN),              ... % PWK ≈ D_n / C_n                   (-)
            "Qn",                       complex(NaN, NaN),              ... % Flow                              (m3/s)
            "Gn",                       complex(NaN, NaN),              ... % Gradient                          (Pa/m)
            "Kn",                       complex(NaN, NaN),              ... % Complex Flow Gain                 (-)
            "tau_n",                    complex(NaN, NaN),              ... % Shear                             (Pa)
            "AnA0",                     complex(NaN, NaN),              ... % Area Puls.                        (-) (m2?)
            "nu_app",                   NaN,                            ... % Viscosity (kinetic)               (m2/s) ? (Pa.s)
            "mu_app",                   NaN,                            ... % Viscosity (dynamic)               (m2/s) ? (Pa.s)
            "H_GQ_n",                   complex(NaN, NaN),              ...
            "H_GQ_n_Geonorm_abs",       NaN,                            ...
            "H_tauQ_n",                 complex(NaN, NaN),              ...
            "H_tauQ_n_Geonorm_abs",     NaN,                            ...
            "H_RQ_n",                   complex(NaN, NaN),              ...
            "H_RQ_n_Geonorm_abs",       NaN,                            ...
            "R_seg_n",                  NaN,                            ...
            "C_seg_n",                  NaN,                            ...
            "P_seg_diss_n",             NaN,                            ...
            "P_seg_store_n",            NaN                             ...
        )                                                           ... 
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
        warning_s("[WOMERSLEY] WomersleyNumberEstimation_n (%i, %i, %i): Not enough valid points in the velocity profile. Skipping fit.", circleIdx, branchIdx, n_harmonic);
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


    D_meas_mag_profile = []; 
    
    % Check if d_profile is provided and valid ([Points x 2 x Time])
    if ~isempty(d_profile) && ndims(d_profile) == 3 && size(d_profile, 3) == numFrames
        
        % We need to process the sliced displacement to match v_profile spatial grid
        % Assuming d_profile was extracted with same spatial points.
        % If v_profile was reduced by valid_idxs, we must reduce d_profile too.
        
        d_profile_reduced = d_profile(valid_idxs, :, :); % Cut background
        
        % Interpolate to match NUM_INTERP_POINTS
        if crossSectionLength > 1
           d_interp = zeros(NUM_INTERP_POINTS, 2, numFrames);
           orig_x = linspace(1, crossSectionLength, crossSectionLength);
           new_x  = linspace(1, crossSectionLength, NUM_INTERP_POINTS);
           
           for t = 1:numFrames
               d_interp(:, 1, t) = interp1(orig_x, squeeze(d_profile_reduced(:, 1, t)), new_x);
               d_interp(:, 2, t) = interp1(orig_x, squeeze(d_profile_reduced(:, 2, t)), new_x);
           end
        else
           d_interp = [];
        end
        
        if ~isempty(d_interp)
             % FFT along time dimension (3)
             d_ft = fftshift(fft(d_interp, N_fft, 3), 3);
             
             % Extract Harmonic for X and Y components separately
             Dx_meas = extractHarmonicProfile(squeeze(d_ft(:, 1, :)), f, cardiac_frequency, n_harmonic);
             Dy_meas = extractHarmonicProfile(squeeze(d_ft(:, 2, :)), f, cardiac_frequency, n_harmonic);
             
             % Compute Magnitude Profile |D_n| = sqrt(|Dx|^2 + |Dy|^2)
             D_meas_mag_profile = sqrt(abs(Dx_meas).^2 + abs(Dy_meas).^2);
        end
    end

    % --- 3. Compute Measured Interaction Statistics (NEW) ---
    if ~isempty(D_meas_mag_profile)
        epsilon = 1e-9;
        
        % Filter based on Velocity signal to avoid noise outside vessel
        abs_v = abs(v_meas);
        abs_D = D_meas_mag_profile;
        
        valid_mask = abs_v > (0.1 * max(abs_v));
        
        if any(valid_mask)
            % 1. Statistic of abs(D_1) (Mean Magnitude in the vessel)
            mean_D_val = mean(abs_D(valid_mask));
            fitParams.metrics.Mean_D1_amp = mean_D_val;
            
            % 2. Ratio abs(D_1) / mean(abs(v_1))
            mean_v_val = mean(abs_v(valid_mask));
            
            fitParams.metrics.Ratio_D1_V1 = mean_D_val / (mean_v_val + epsilon);
        end
    end

    x_coords = linspace(-1, 1, NUM_INTERP_POINTS);
    
    FIXED_NU = 3e-6; % m^2/s (Phase 1 recommendation)
    OMEGA_N = 2 * pi * cardiac_frequency * n_harmonic;
    
    % Factor to convert 'width' (normalized 0-1) to Radius in Meters
    % PIXEL_SIZE is in microns, convert to meters (1e-6)
    % R_SCALE_METERS = (PIXEL_SIZE * 1e-6) * crossSectionLength / 2; 
    
    % Setup all parameters for fit depending on options
    % p = [real(Cn), imag(Cn), real(Dn), imag(Dn), center]
    Cn_init_complex = mean(v_meas(abs(v_meas) > 0));
    Dn_init_complex = 0; % Start with no wall motion
    center_init     = init_fit.geoParams.center_norm;
    [p_init, lb, ub, modelArgHandle] = fit_setup(options, Cn_init_complex, Dn_init_complex, center_init, FIXED_NU);

    % Which is basically init_fit.geoParams.width_fit, the fit on R0
    width_norm = fitParams.width;

    % costFunctionHandle = @(p) costFun(p, x_coords, v_meas, width_norm, PSF_KERNEL, FIXED_NU, R_SCALE_METERS, OMEGA_N, fitParams.R0);
    costFunctionHandle = @(p) costFun(p, modelArgHandle, x_coords, v_meas,  width_norm, PSF_KERNEL, fitParams.R0, OMEGA_N);

    fit_options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective');
    try
        [p_fit, ~, residual, exitflag, ~, ~, jacobian] = lsqnonlin(costFunctionHandle, p_init, lb, ub, fit_options);
        [Cn_fit, Dn_fit, center_fit, nu_fit] = modelArgHandle(p_fit);

        ALPHA_N = fitParams.R0 * sqrt(OMEGA_N / nu_fit);
        %          generate_womersley_profile(x, Cn, Dn, center, alpha, width_norm, psf_kernel)
        uWom_fit = generate_womersley_profile(x_coords, Cn_fit, Dn_fit, center_fit, ALPHA_N, width_norm, PSF_KERNEL);

        [Kappa_n, res_mag_RMS, res_phase_RMS, res_phase_RMS_msk, harmonic_SNR_dB] = calculate_fit_precision_metrics(v_meas, residual, jacobian);

        fitParams.Kappa_n                = Kappa_n;
        fitParams.residual_mag_RMS       = res_mag_RMS;    
        fitParams.residual_phase_RMS     = res_phase_RMS;
        fitParams.residual_phase_RMS_msk = res_phase_RMS_msk;
        fitParams.harmonic_SNR_dB        = harmonic_SNR_dB;
        fitParams.fit_exitflag           = exitflag;
        
        % fitParams.Kappa_n = sqrt(condest(jacobian'*jacobian));

        % Cn_fit              = p_fit(1) + 1i * p_fit(2);
        % Dn_fit              = p_fit(3) + 1i * p_fit(4);
        % center_fit          = p_fit(5);
    
        fitParams.alpha_n   = ALPHA_N;
        fitParams.Cn        = Cn_fit;
        fitParams.Dn        = Dn_fit;
        fitParams.center    = center_fit;

        % Already done above, but for clarity
        fitParams.omega_0 = 2 * pi * cardiac_frequency;
        fitParams.omega_n = fitParams.omega_0 * n_harmonic;

        fitParams.harmonic = n_harmonic;
    
        if options.ModelType == "rigid"
            fitParams.Rn = complex(NaN, NaN);
        else
            fitParams.Rn = (fitParams.Dn / fitParams.Cn) * fitParams.R0;
        end
    
        fitParams.metrics.nu_app = nu_fit;

        fitParams.metrics = calculate_symbols(fitParams, RHO_BLOOD, options, D_meas_reg);

    catch ME 
        warning_s("[WOMERSLEY] Womersley fit failed for %s (idx %d): %s", name, idx, ME.message);
        return;
    end
    
    % [parabole_fit_systole, parabole_fit_diastole] = analyse_lumen_size(v_profile_good_idx_sav, SYS_IDXS, DIAS_IDXS);
    % 
    % estimated_width.systole = parabole_fit_systole;
    % estimated_width.diastole = parabole_fit_diastole;
    
    % ============================ [ Figures ] ========================== %
    if params.saveFigures
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
        line2_str = sprintf('Center: %-10.2f     |R_n/R_0| : %.2f %%', fitParams.center, abs(fitParams.metrics.RnR0_complex * 100));
        line3_str = sprintf('Width : %-10.2f     Phase(R_n): %.1f°', fitParams.width, rad2deg(angle(fitParams.metrics.RnR0_complex)));
        
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
            warning("[WOMERSLEY] Could not save figure");
        end
        
        if ~strcmpi(get(hFig, 'Visible'), 'on')
            close(hFig);
        end
    end
end


% +=====================================================================+ %
% |                          HELPER FUNCTIONS                           | %
% +=====================================================================+ %

% ============================ [ FIT SETUP ] ============================ %

function [p_init, lb, ub, modelArgHandle] = fit_setup(options, Cn_init, Dn_init, center_init, FIXED_NU)
    arguments
        options, Cn_init, Dn_init, center_init, FIXED_NU
    end

    if options.ModelType == "rigid" && options.NuType == "fixed"
        p_init = [real(Cn_init), imag(Cn_init), center_init];
        
        lb = [-Inf, -Inf, -0.8];
        ub = [ Inf,  Inf,  0.8];

        modelArgHandle = @(p) deal(complex(p(1), p(2)), complex(0, 0), p(3), FIXED_NU);

    elseif options.ModelType == "moving" && options.NuType == "fixed"
        p_init = [real(Cn_init), imag(Cn_init), real(Dn_init), imag(Dn_init), center_init];
        
        lb = [-Inf, -Inf, -Inf, -Inf, -0.8];
        ub = [ Inf,  Inf,  Inf,  Inf,  0.8];

        modelArgHandle = @(p) deal(complex(p(1), p(2)), complex(p(3), p(4)), p(5), FIXED_NU);
    elseif options.ModelType == "rigid" && options.NuType == "free"
        p_init = [real(Cn_init), imag(Cn_init), center_init, FIXED_NU];
        
        lb = [-Inf, -Inf, -0.8, 1.5e-6];
        ub = [ Inf,  Inf,  0.8, 4.0e-6];

        modelArgHandle = @(p) deal(complex(p(1), p(2)), complex(0, 0), p(3), p(4));
    else
        warning_s("[WOMERSLEY] fit_setup: Invalid options combinaison (%s, %s)", options.ModelType, options.NuType);
    end
end

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
         warning("[WOMERSLEY] fitGeometryFromMean: Profile too short"); 
         return; 
    end

    % TODO: Maybe not fit on interpolated values
    v_profile_interp = interp1(linspace(-1, 1, crossSectionLength), v_profile, x_grid_normalized);
    
    % PSF Kernel
    PIXEL_SIZE = params.px_size * 1e-3;

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

function res = costFun(p, modelArgHandle, x_coords, v_meas, width_norm, psf_kernel, R0, omega_n)
    [Cn, Dn, center, nu] = modelArgHandle(p);

    % Based on the equation 18a of paper v7
    alpha_n = R0 * sqrt(omega_n / nu);

    val = generate_womersley_profile(x_coords, Cn, Dn, center, alpha_n, width_norm, psf_kernel);

    res = [real(val - v_meas.'); ...
           imag(val - v_meas.')];
end

function Bn = uWom_base(alpha, r)
    lambda_n = 1i^(3/2) * alpha;
    J0 = besselj(0, lambda_n);

    if J0 == 0
        Bn = zeros(size(r));
        return;
    end

    Bn = (1 - (besselj(0, lambda_n * r) ./ J0));
end


function Psin = uWom_psi(alpha, r)
    % Formula: Psi_n = - [lambda * J1(lambda) / J0(lambda)^2] * J0(lambda*r)
    %
    % Inputs:
    %   alpha - The Womersley number (a dimensionless real number).
    %   r     - The normalized radial coordinate (r/R0), typically a vector 
    %           from 0 to 1.
    
    % it is ((-1 + i) / sqrt(2)) * alpha
    lambda_n = 1i^(3/2) * alpha;

    numerator_scalar = -lambda_n * besselj(1, lambda_n);
    denominator_scalar = besselj(0, lambda_n) ^ 2;

    if denominator_scalar == 0
        Psin = zeros(size(r));
        return;
    end

    complex_scalar = numerator_scalar / denominator_scalar;
    
    radial_profile = besselj(0, lambda_n .* r);
    
    Psin = complex_scalar .* radial_profile;
end

% ============================= [ WOMERSELY ] =========================== %

function model_profile = generate_womersley_profile(x, Cn, Dn, center, alpha, width_norm, psf_kernel)
    % A more generalized and clean version

    arguments
        x, Cn, Dn, center, alpha, width_norm, psf_kernel
    end

    % Normalise coordinates
    r = (x - center) / width_norm;

    [Bn, Psin] = get_womersley_basis(alpha, r);

    if ~isempty(psf_kernel)
        Bn = conv(Bn, psf_kernel, 'same');
        if Dn ~= 0
            Psin = conv(Psin, psf_kernel, 'same');
        end
    end

    model_profile = (Cn * Bn) + (Dn * Psin);
    
    model_profile(abs(r) > 1) = 0;
end


function [Bn, Psin] = get_womersley_basis(alpha, r)
    arguments
        alpha
        r
    end
    % Returns the Womersley basis as calculated
    % C.f.  18a - 18c of v7

    Bn = uWom_base(alpha, r);

    if nargout > 1
        Psin = uWom_psi(alpha, r);
    else
        Psin = [];
    end
end

function model_profile = generate_moving_wall_model(p, x, width_norm, psf_kernel, nu, r_scale, omega_n, R0)
    % psf_kernel is optionnal

    % p = [alpha, real(Cn), imag(Cn), real(Dn), imag(Dn), center, width]
    % alpha_1 = p(1);
    Cn      = p(1) + 1i * p(2);
    Dn      = p(3) + 1i * p(4);
    center  = p(5);

    width   = width_norm;
    
    r = (x - center) / width;


    % See paper section 4.1, the alpha can be deduced from alpha_1
    % The functions Bn abd Psin use an alpha that is dimensionless
    % (calculated, before hand)
    % alpha_n = alpha_1 * sqrt(n_harmonic);
    
    % R = width * r_scale; % (Rn)
    R = R0; % (R0) 

    if R <= 0
        alpha_n = 1e-3; 
    else
        alpha_n = R * sqrt(omega_n / nu);
    end
    
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
    
    model_profile = amplitude * profile;
end

% ======================================================================= %
    
function v_meas = extractHarmonicProfile(v_profile_ft, f_vector, base_frequency, n_harmonic, dimToMean)
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

    arguments
        v_profile_ft, f_vector, base_frequency, n_harmonic
        dimToMean = 2
    end

    target_frequency = n_harmonic * base_frequency;

    freq_resolution = f_vector(2) - f_vector(1);
    margin_hz = freq_resolution * 1.5;
    
    harmonic_indices = find(abs(f_vector - target_frequency) <= margin_hz);

    if isempty(harmonic_indices)
        warning("[WOMERSLEY] Harmonic %d (%.2f Hz) not found in frequency band. Using closest single peak instead.", n_harmonic, target_frequency);
        [~, closest_idx] = min(abs(f_vector - target_frequency));
        harmonic_indices = closest_idx;
    end

    v_meas = mean(v_profile_ft(:, harmonic_indices), dimToMean);
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
        warning("[WOMERSLEY] Cannot create PSF kernel: The spatial resolution (dx_um) is zero.");
    end
    
    fwhm_in_points = fwhm_um / dx_um;

    sigma_in_points = fwhm_in_points / (2 * sqrt(2 * log(2)));

    kernel_radius = ceil(3 * sigma_in_points);
    kernel_size = 2 * kernel_radius + 1;

    % 'fspecial' creates a normalized kernel that sums to 1
    psf_kernel = fspecial('gaussian', [1, kernel_size], sigma_in_points);
end


% ============================== [ SYMBOLS ] ============================ %

function [Kappa_n, res_mag_RMS, res_phase_RMS, res_phase_RMS_msk, harmonic_SNR_dB] = calculate_fit_precision_metrics(v_meas, residual, jacobian)

    v_meas   = v_meas(:);
    residual = residual(:);

    % ================ KAPPA_N ================
    Kappa_n = cond(full(jacobian)); % sqrt(condest(jacobian'*jacobian));

    % ================   RMS   ================

    N_points = numel(v_meas);

    % Résidus renvoyés par lsqnonlin : [Re(eps); Im(eps)]
    eps_complex = residual(1:N_points) + 1i * residual(N_points+1:2*N_points);

    % Données et modèle au point optimum
    v_data_col  = v_meas(:);
    v_model_col = v_data_col - eps_complex;

    % Résidu amplitude RMS normalisé
    mag_data  = abs(v_data_col);
    mag_model = abs(v_model_col);

    diff_mag       = mag_data - mag_model;
    diff_mag_sq    = diff_mag .^ 2;
    mse_mag        = mean(diff_mag_sq);       % Mean Squared Error
    rms_data       = sqrt(mean(mag_data.^2)); % RMS du signal

    res_mag_RMS = sqrt(mse_mag) / (rms_data + eps);  % eviter div/0

    % Résidu phase RMS
    % phase_err = angle(V_data * conj(V_model)) dans [-pi, pi]
    phase_diff_rad = angle(v_data_col .* conj(v_model_col));

    % ponderation par amplitude pour ignorer le bruit
    mask = mag_data > 0.1 * max(mag_data);
    res_phase_RMS_msk = sqrt(mean(phase_diff_rad(mask) .^ 2));
    res_phase_RMS     = sqrt(mean(phase_diff_rad       .^ 2));

    % ================   SNR   ================

    signal_power = sum(abs(v_model_col) .^ 2);
    noise_power  = sum(abs(eps_complex) .^ 2);

    harmonic_SNR_dB = 10 * log10(max(signal_power, eps) / max(noise_power, eps));
end


function K = calculate_K_factor(lambda)
    % Calculates the flow moment factor K(alpha).
    % K(alpha) = 1 - (2*J1(lambda) / (lambda * J0(lambda)))
    
    if lambda == 0 || besselj(0, lambda) == 0
        K = complex(NaN, NaN); % Avoid division by zero
        return;
    end
    
    K = 1 - (2 * besselj(1, lambda)) / (lambda * besselj(0, lambda));
end


function metrics = calculate_symbols(fitParams, rho_blood, options, D_meas_reg)
    Cn          = fitParams.Cn;
    Dn          = fitParams.Dn;
    R0          = fitParams.R0; % in meters
    % alpha_1     = fitParams.alpha_1;
    alpha_n     = fitParams.alpha_n;
    omega_n     = fitParams.omega_n;

    Rn          = fitParams.Rn;

    lambda_n = (1i^(3/2)) * alpha_n;

    metrics = struct();
    % metrics = struct(...
    %     'RnR0_complex',     NaN, ... % PWK ≈ D_n / C_n
    %     'RnR0_mag',         NaN, ...
    %     'RnR0_phase_deg',   NaN, ...
    %     'Qn',              NaN, ... % Flow (mm3/s)
    %     'Gn',              NaN, ... % Gradient (Pa/m)
    %     'tau_n',            NaN, ... % Shear (Pa)
    %     'H_Zn',             NaN, ... % Impedance (Pa s/mm3/m)
    %     'H_Rn',             NaN, ... % Geom-Norm (Pa)
    %     'AnA0',             NaN, ...
    %     'nu_app',           NaN, ... % Viscosity (kinetic) (m2/s)
    %     'mu_app',           NaN  ... % Viscosity (dinamic)
    % );

    % Flow (Qn)
    metrics.Kn = calculate_K_factor(lambda_n);
    metrics.Qn = pi * R0 ^ 2 * Cn * metrics.Kn; % units: m³/s

    % Gradient (Gn)
    metrics.Gn = 1i * omega_n * rho_blood * Cn; % units: Pa/m

    % Per length impedance ratio (Pa.s/m4)
    metrics.H_GQ_n = metrics.Gn / metrics.Qn;
    % Geometry normalized impedance magnitude (Pa.s/m2)
    metrics.H_GQ_n_Geonorm_abs = norm(metrics.H_GQ_n) * (R0 ^ 2);


    % Handle Moving Wall metrics to be not calculated when Rigid
    if options.ModelType == "rigid"
        metrics.RnR0_complex        = complex(NaN, NaN);
        metrics.AnA0                = complex(NaN, NaN);
        metrics.H_RQ_n              = complex(NaN, NaN);
        metrics.H_RQ_n_Geonorm_abs  = NaN;
    else
        % PWK (Rn/R0)
        if Cn ~= 0
            metrics.RnR0_complex = Dn / Cn;
        else
            metrics.RnR0_complex = complex(NaN, NaN);
        end

        % Area puls. (An/A0)
        metrics.AnA0 = 2 * metrics.RnR0_complex;

        % Dilatation per unit pulsatile flow (s/m3)
        metrics.H_RQ_n = (Rn / R0) / metrics.Qn;
        metrics.H_RQ_n_Geonorm_abs = norm(metrics.H_RQ_n) * (R0 ^ 2);
    end

    % Viscosity (νapp)
    % numerator = (R0 ^ 2) * fitParams.omega_0;
    % denominator = alpha_1 ^ 2;
    
    % if denominator > 0
    %     metrics.nu_app = numerator / denominator;
    %     metrics.mu_app = metrics.nu_app * rho_blood;
    % else
    %     metrics.nu_app = NaN;
    %     metrics.mu_app = NaN;
    % end

    metrics.nu_app = fitParams.metrics.nu_app;
    metrics.mu_app = metrics.nu_app * rho_blood;
    
    % Shear (τn)
    if ~isnan(metrics.mu_app) && R0 > 0 && besselj(0, lambda_n) ~= 0
        metrics.tau_n = (metrics.mu_app * Cn * lambda_n * besselj(1, lambda_n)) / (R0 * besselj(0, lambda_n)); % units: Pa
    else
        metrics.tau_n = complex(NaN, NaN);
    end

    % Shear stress per unit pulsatile flow (Pa.s/m3)
    metrics.H_tauQ_n = metrics.tau_n / metrics.Qn;
    metrics.H_tauQ_n_Geonorm_abs = norm(metrics.H_tauQ_n) * (R0 ^ 2);

    % Resistance-like contribution per unit length (Pa.s/m4)
    metrics.R_seg_n = real(metrics.H_GQ_n);

    % Compliance-like (inverse-compliance) per unit length (Pa.s/m4)
    metrics.C_seg_n = 1 / imag(metrics.H_GQ_n);

    % Dissipative AC power per unit length (W/m)
    metrics.P_seg_diss_n  = 1/2 * (norm(metrics.Qn) ^ 2) * real(metrics.H_GQ_n); 

    % Reactive (stored) AC power per unit length (W/m)
    metrics.P_seg_store_n = 1/2 * (norm(metrics.Qn) ^ 2) * imag(metrics.H_GQ_n);

end

function h_metrics = calculate_harmonic_metrics(fitParams)
    arguments
        fitParams % list of fitParams based uppon the harmonics
    end

    fitParams = fitParams(:);
    leng_harmo = size(fitParams, 1);

    metrics = [fitParams.metrics];

    h_metrics = struct();

    % Shear harmonic ratios |τ2|/|τ1|, |τ3|/|τ1| (-)
    h_metrics.RhoTau21 = ternary_op(leng_harmo >= 2, @() norm(metrics(2).tau_n) / norm(metrics(1).tau_n), @() NaN);
    h_metrics.RhoTau31 = ternary_op(leng_harmo >= 3, @() norm(metrics(3).tau_n) / norm(metrics(1).tau_n), @() NaN);

    % Phase skewness of shear, ϕτ,2 − 2ϕτ,1 (°)
    h_metrics.DeltaPhiTau2 = ternary_op(leng_harmo >= 2, @() angle(metrics(2).tau_n) - 2 * angle(metrics(1).tau_n), @() NaN);

    % TODO: Add Description, because I can't Find anything
    h_metrics.RhoQ21       = ternary_op(leng_harmo >= 2, @() norm(metrics(2).Qn) / norm(metrics(1).Qn), @() NaN);
    h_metrics.RhoQ31       = ternary_op(leng_harmo >= 3, @() norm(metrics(3).Qn) / norm(metrics(1).Qn), @() NaN);
    h_metrics.DeltaPhiQ2   = ternary_op(leng_harmo >= 2, @() angle(metrics(2).Qn) - 2 * angle(metrics(1).Qn), @() NaN);
end

function derived = calculate_segment_derived_metrics(fitParams, qc_model, geoParams, rho_blood)
% CALCULATE_SEGMENT_DERIVED_METRICS
%   Compute the derived metrics defined in "Derived Womersley Metrics for
%   Retinal Doppler Holography" for a single vessel segment and one model.
%
%   INPUT
%     fitParams  : [nH x 1] struct array (same as in segments_metrics.<Model>)
%     qc_model   : QC struct returned by womersley_qc for this model
%                  (fields: harmonic_valid, use_ratios, reason)
%     geoParams  : geometry struct from fitGeometryOnMean (R0, DC_Amp, ...)
%     rho_blood  : blood density (kg/m^3)
%
%   OUTPUT
%     derived    : struct with fields:
%                  Q0, tau0_mag, PI_Q, PI_V, PI_tau,
%                  phi_GQ1, phi_tauQ1,
%                  LossFactor1, ReactanceFactor1,
%                  PowerStoreOverDiss1, PowerStoreFraction1,
%                  HF_flow_index, HF_shear_index
    arguments
        fitParams
        qc_model
        geoParams
        rho_blood
    end

    nH = numel(fitParams);

    derived = struct();
    derived.Q0                  = NaN;
    derived.tau0_mag            = NaN;
    derived.PI_Q                = NaN;
    derived.PI_V                = NaN;
    derived.PI_tau              = NaN;
    derived.phi_GQ1             = NaN;
    derived.phi_tauQ1           = NaN;
    derived.LossFactor1         = NaN;
    derived.ReactanceFactor1    = NaN;
    derived.PowerStoreOverDiss1 = NaN;
    derived.PowerStoreFraction1 = NaN;
    derived.HF_flow_index       = NaN;
    derived.HF_shear_index      = NaN;

    if nH == 0
        return;
    end

    % ==================== DC quantities: Q0 and |tau0| ===================== %

    R0 = geoParams.R0;              % (m)
    V0_center = geoParams.DC_Amp;   % DC centreline amplitude (complex)

    if isnan(R0) || R0 <= 0 || isnan(V0_center)
        Q0       = NaN;
        Q0_mag   = NaN;
        tau0_mag = NaN;
    else
        % v0(r) = V0_center * (1 - (r^2 / R0^2))
        % V0_mean = V0_center / 2
        % Q0 = pi * R0^2 * V0_mean = (pi/2) * R0^2 * V0_center
        Q0     =  (pi * 0.5) * R0^2 * V0_center;
        Q0_mag = abs(Q0);

        % choose nu_app from first valid harmonic
        nu_app = NaN;
        valid_idx = find(qc_model.harmonic_valid(:));

        if ~isempty(valid_idx)
            % try to pick harmonic 1 if available
            allH = [fitParams.harmonic];
            % Really hugly, but searches for harmonics == 1 among valid_idx
            idx1 = valid_idx(find(allH(valid_idx) == 1, 1, "first"));
            if isempty(idx1)
                idx1 = valid_idx(1);
            end
            nu_app = fitParams(idx1).metrics.nu_app;
        end

        if isnan(nu_app) || nu_app <= 0
            nu_app = 3e-6;  % fallback
        end

        mu_app   = rho_blood * nu_app;
        V0c_mag  = abs(V0_center);

        % |tau0| = mu_app * 2|V0_center| / R0
        if V0c_mag > 0
            tau0_mag = mu_app * (2 * V0c_mag / R0);
        else
            tau0_mag = NaN;
        end
    end

    derived.Q0       = Q0;
    derived.tau0_mag = tau0_mag;

    % =============== Grab first harmonic quantities (n = 1) ================ %

    allH = [fitParams.harmonic];
    idx1 = find(allH == 1, 1, "first");
    if isempty(idx1) || ~qc_model.harmonic_valid(idx1)
        % No valid first harmonic
        derived.PI_Q                = NaN;
        derived.PI_V                = NaN;
        derived.PI_tau              = NaN;
        derived.phi_GQ1             = NaN;
        derived.phi_tauQ1           = NaN;
        derived.LossFactor1         = NaN;
        derived.ReactanceFactor1    = NaN;
        derived.PowerStoreOverDiss1 = NaN;
        derived.PowerStoreFraction1 = NaN;
        derived.HF_flow_index       = NaN;
        derived.HF_shear_index      = NaN;
        return;
    end

    fp1      = fitParams(idx1);
    C1       = fp1.Cn;
    Q1       = fp1.metrics.Qn;
    tau1     = fp1.metrics.tau_n;
    HGQ1     = fp1.metrics.H_GQ_n;
    HtauQ1   = fp1.metrics.H_tauQ_n;
    Pdiss1   = fp1.metrics.P_seg_diss_n;
    Pstore1  = fp1.metrics.P_seg_store_n;

    Q1_mag   = abs(Q1);
    tau1_mag = abs(tau1);

    % =============== Pulsatility indices PI_Q, PI_V, PI_tau ================ %

    if ~isnan(Q0_mag) && Q0_mag > 0 && Q1_mag > 0
        derived.PI_Q = Q1_mag / Q0_mag;
    else
        derived.PI_Q = NaN;
    end

    V0c_mag = abs(V0_center);
    if V0c_mag > 0 && ~isnan(V0c_mag) && abs(C1) > 0
        derived.PI_V = abs(C1) / V0c_mag;
    else
        derived.PI_V = NaN;
    end

    if ~isnan(tau0_mag) && tau0_mag > 0 && tau1_mag > 0
        derived.PI_tau = tau1_mag / tau0_mag;
    else
        derived.PI_tau = NaN;
    end

    % =========== Impedance phases and factors (Loss / Reactance) =========== %

    if ~isnan(HGQ1) && HGQ1 ~= 0
        HGQ1_mag = abs(HGQ1);
        derived.phi_GQ1          = angle(HGQ1);
        derived.LossFactor1      = real(HGQ1) / HGQ1_mag;
        derived.ReactanceFactor1 = imag(HGQ1) / HGQ1_mag;
    else
        derived.phi_GQ1          = NaN;
        derived.LossFactor1      = NaN;
        derived.ReactanceFactor1 = NaN;
    end

    if ~isnan(HtauQ1) && HtauQ1 ~= 0
        derived.phi_tauQ1 = angle(HtauQ1);
    else
        derived.phi_tauQ1 = NaN;
    end

    % ======================== Power balance metrics ======================== %

    if ~isnan(Pdiss1) && Pdiss1 ~= 0
        derived.PowerStoreOverDiss1 = Pstore1 / Pdiss1;
    else
        derived.PowerStoreOverDiss1 = NaN;
    end

    if ~isnan(Pdiss1 + Pstore1) && (Pdiss1 + Pstore1) ~= 0
        derived.PowerStoreFraction1 = Pstore1 / (Pstore1 + Pdiss1);
    else
        derived.PowerStoreFraction1 = NaN;
    end

    % ======= Higher harmonic content indices (HF_flow and HF_shear) ======== %
    % Only include harmonics that pass QC and have n >= 2

    if ~qc_model.use_ratios || Q1_mag <= 0 || tau1_mag <= 0
        derived.HF_flow_index  = NaN;
        derived.HF_shear_index = NaN;
        return;
    end

    valid_hf = qc_model.harmonic_valid(:) & (allH(:) >= 2);

    if ~any(valid_hf)
        derived.HF_flow_index  = 0;   % no higher harmonics
        derived.HF_shear_index = 0;
        return;
    end

    metrics_all = [fitParams.metrics];
    Qn_all      = [metrics_all.Qn];
    taun_all    = [metrics_all.tau_n];

    Qn_hf   = Qn_all(valid_hf);
    taun_hf = taun_all(valid_hf);

    num_flow  = sum(abs(Qn_hf).^2);
    num_shear = sum(abs(taun_hf).^2);

    derived.HF_flow_index  = num_flow  / (Q1_mag^2);
    derived.HF_shear_index = num_shear / (tau1_mag^2);
end


% ================================ [ QC ] =============================== %

function qc = womersley_qc(fitParams, modelType, varargin)
% WOMERSLEY_QC  Simple QC gating around a fitParams struct array
%
%   qc = womersley_qc(fitParams, modelType)
%
%   INPUT
%     fitParams : [nHarmonics x 1] struct array from WomersleyNumberEstimation_n
%     modelType : "rigid" or "moving" (same as options.ModelType)
%
%   OPTIONAL name/value thresholds:
%     'SNR_min'          : minimum harmonic_SNR_dB    (default 5 dB)
%     'log10Kappa_max'   : max log10(Kappa_n)         (default: 4 (rigid), 6 (moving))
%     'max_RnR0'         : max |Rn/R0| for moving wall (default 0.3)
%     'min_Q1_rel'       : min |Q1| relative to max |Qn| for ratios (default 0.1)
%
%   OUTPUT
%     qc.harmonic_valid  : [nH x 1] logical valid/invalid per harmonic
%     qc.use_ratios      : scalar logical (true if it is safe to use RhoQ21, RhoTau21, ...)
%     qc.reason          : struct of masks for debugging (fields: snr, kappa, wall, success, q1)

    arguments
        fitParams
        modelType   string {mustBeMember(modelType, ["rigid","moving"])}
    end
    arguments (Repeating)
        varargin
    end

    % --- Default thresholds ---
    SNR_min        = 5;    % dB
    if modelType == "moving"
        log10Kappa_max = 6;    % 1e6
    else
        log10Kappa_max = 4;    % 1e4
    end
    max_RnR0       = 0.3;  % 30% radius pulsation
    min_Q1_rel     = 0.1;  % |Q1| >= 10% of max |Qn| in this segment

    % --- Override via name/value pairs if needed ---
    if ~isempty(varargin)
        for i = 1:2:numel(varargin)
            name  = varargin{i};
            value = varargin{i+1};
            switch lower(name)
                case 'snr_min'
                    SNR_min = value;
                case 'log10kappa_max'
                    log10Kappa_max = value;
                case 'max_rnr0'
                    max_RnR0 = value;
                case 'min_q1_rel'
                    min_Q1_rel = value;
                otherwise
                    warning('womersley_qc: unknown parameter "%s"', name);
            end
        end
    end

    nH = numel(fitParams);
    if nH == 0
        qc.harmonic_valid = false(0,1);
        qc.use_ratios     = false;
        qc.reason         = struct();
        return;
    end

    % --- Basic per-harmonic masks ---
    snr_ok     = false(nH,1);
    kappa_ok   = false(nH,1);
    success_ok = false(nH,1);
    wall_ok    = true(nH,1);   % only used for moving wall

    for k = 1:nH
        fp = fitParams(k);

        % Success of solver
        success_ok(k) = (isfield(fp, "fit_exitflag") && fp.fit_exitflag > 0);

        % SNR gate
        if isfield(fp, "harmonic_SNR_dB") && ~isnan(fp.harmonic_SNR_dB)
            snr_ok(k) = (fp.harmonic_SNR_dB >= SNR_min);
        end

        % Conditioning gate
        if isfield(fp, "Kappa_n") && ~isnan(fp.Kappa_n) && fp.Kappa_n > 0
            kappa_ok(k) = (log10(fp.Kappa_n) <= log10Kappa_max);
        else
            kappa_ok(k) = false;
        end

        % Moving-wall: limit |Rn/R0|
        if modelType == "moving"
            if isfield(fp, "metrics") && isfield(fp.metrics, "RnR0_complex")
                val = fp.metrics.RnR0_complex;
                if ~isnan(val)
                    wall_ok(k) = (abs(val) <= max_RnR0);
                else
                    wall_ok(k) = false;
                end
            else
                wall_ok(k) = false;
            end
        end
    end

    harmonic_valid = snr_ok & kappa_ok & success_ok & wall_ok;

    % --- Ratio safety: is it safe to trust |Q2|/|Q1|, |τ2|/|τ1|, etc.? ---
    % Require:
    %   - harmonic 1 is valid
    %   - |Q1| not "tiny" compared to max |Qn| in this segment
    Qn = arrayfun(@(fp) fp.metrics.Qn, fitParams);
    Qabs = abs(Qn);

    if nH >= 1 && harmonic_valid(1) && any(Qabs > 0)
        Q1   = Qabs(1);
        Qmax = max(Qabs);
        q1_ok = (Q1 >= min_Q1_rel * Qmax);
    else
        q1_ok = false;
    end

    qc.harmonic_valid = harmonic_valid;
    qc.use_ratios     = q1_ok;

    % Optional debug info
    qc.reason = struct( ...
        "snr",     snr_ok, ...
        "kappa",   kappa_ok, ...
        "success", success_ok, ...
        "wall",    wall_ok, ...
        "q1_ok",   q1_ok ...
    );
end


% +=====================================================================+ %
% |                                DEBUG                                | %
% +=====================================================================+ %

function visualizeDCFit(v_mean, geoParams, psf_kernel) %#ok<DEFNU>
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

  