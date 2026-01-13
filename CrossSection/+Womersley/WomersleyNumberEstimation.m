function results = WomersleyNumberEstimation(v_profile, cardiac_frequency, name, circleIdx, branchIdx, d_profile)
    arguments
        v_profile
        cardiac_frequency
        name
        circleIdx
        branchIdx
        d_profile
    end

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
    FIXED_NU = 3e-6; % m^2/s (Phase 1 recommendation)

    psf_kernel = createGaussianPSFKernel(FWHM_um, NUM_INTERP_POINTS, crossSectionLength, PIXEL_SIZE);
    % Fit a simple PSF-convolved Parabolic/Plug model to get Geometry
    [geoParams, ~] = fitGeometryOnMean(v_profile, psf_kernel, ToolBox);
    DC_metrics = calculateDCmetrics(geoParams, RHO_BLOOD, FIXED_NU);

    if isnan(geoParams.R0)
        warning("[WOMERSLEY] Geometry fit failed");
        return;
    end

    init_fit = struct(...
        "psf_kernel",   psf_kernel, ...
        "geoParams",    geoParams   ...
    );

    HARMONIC_NUMBER = params.json.exportCrossSectionResults.Womersley.MaxHarmonic;

    % TODO: Parfor does not seem to work with toolbox
    for i = 1:HARMONIC_NUMBER
        fitParams.RigidWallFixedNu(i)  = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, circleIdx, branchIdx, i, init_fit, d_profile, FIXED_NU, ToolBox, ModelType="rigid");
        fitParams.MovingWallFixedNu(i) = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, circleIdx, branchIdx, i, init_fit, d_profile, FIXED_NU, ToolBox, ModelType="moving");
        fitParams.RigidWallFreeNu(i)   = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, circleIdx, branchIdx, i, init_fit, d_profile, FIXED_NU, ToolBox, ModelType="rigid", NuType="free");
    end

    % Quality Control
    qc.qc_RigidWallFixedNu  = calculateWomersleyQC(fitParams.RigidWallFixedNu,  "rigid");
    qc.qc_MovingWallFixedNu = calculateWomersleyQC(fitParams.MovingWallFixedNu, "moving");
    qc.qc_RigidWallFreeNu   = calculateWomersleyQC(fitParams.RigidWallFreeNu,   "rigid");

    % h_metrics.RigidWallFixedNu  = calculateHarmonicMetrics(fitParams.RigidWallFixedNu);
    % h_metrics.MovingWallFixedNu = calculateHarmonicMetrics(fitParams.MovingWallFixedNu);
    % h_metrics.RigidWallFreeNu   = calculateHarmonicMetrics(fitParams.RigidWallFreeNu);

    % filtered metrics (Should be done inside Analysis Program, not Here)
    h_metrics.RigidWallFixedNu  = calculateHarmonicMetrics(fitParams.RigidWallFixedNu(qc.qc_RigidWallFixedNu.harmonic_valid));
    h_metrics.MovingWallFixedNu = calculateHarmonicMetrics(fitParams.MovingWallFixedNu(qc.qc_MovingWallFixedNu.harmonic_valid));
    h_metrics.RigidWallFreeNu   = calculateHarmonicMetrics(fitParams.RigidWallFreeNu(qc.qc_RigidWallFreeNu.harmonic_valid));

    derived.RigidWallFixedNu  = calculateSegmentDerivedMetrics(fitParams.RigidWallFixedNu, qc.qc_RigidWallFixedNu, init_fit.geoParams, RHO_BLOOD);
    derived.MovingWallFixedNu = calculateSegmentDerivedMetrics(fitParams.MovingWallFixedNu, qc.qc_MovingWallFixedNu, init_fit.geoParams, RHO_BLOOD);
    derived.RigidWallFreeNu   = calculateSegmentDerivedMetrics(fitParams.RigidWallFreeNu, qc.qc_RigidWallFreeNu, init_fit.geoParams, RHO_BLOOD);

    % metrics for each segments / harmonic
    results.segments_metrics = fitParams;

    % metrics for each segments (using multiple harmonics)
    results.harmonic_metrics = h_metrics;

    % QCs are both, mostly calculated for each segments / harmonics, but some 
    % are just for segments
    results.qc = qc;
    results.DC_metrics = DC_metrics;

    % TODO: Maybe removed derived and put inside respective results
    % Derived are metrics that are some extra metrics, they could be included 
    % inside segment_metrics and harmonic
    % Some are segments/harmonics and some not
    results.derived = derived;
end


function fitParams = WomersleyNumberEstimation_n(v_profile, cardiac_frequency, name, circleIdx, branchIdx, n_harmonic, init_fit, d_profile, FIXED_NU, ToolBox, options)
    arguments
        v_profile, cardiac_frequency, name, circleIdx, branchIdx, n_harmonic, init_fit, d_profile, FIXED_NU, ToolBox
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
    
    % Really should be a power of 2
    FFT_PADDING_FACTOR = 16;

    RHO_BLOOD = 1060; % Density of blood in kg/m^3

    PSF_KERNEL = init_fit.psf_kernel;

    fitParams = getFitParamsStruct();
    fitParams.R0    = init_fit.geoParams.R0;
    fitParams.width = init_fit.geoParams.width_norm;
    
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
    % Will pad to the next power of 2 to ease the FFT algorithm
    N_fft = (2 ^ nextpow2(numFrames)) * FFT_PADDING_FACTOR;
    v_profile_ft = fftshift(fft(v_profile, N_fft, 2), 2);
    
    f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, N_fft);
    
    % fitParams.v_profile_ft = v_profile_ft;
    % fitParams.frequency_vector = f;

    % figure;
    % subplot(2,1,1);
    % plot(mean(abs(fft(v_profile, N_fft, 2)), 1)); title("FFT v\_profile with PADDING 16"); xlim([0, N_fft]);
    % subplot(2,1,2);
    % plot(mean(abs(fft(v_profile, (2 ^ nextpow2(numFrames)), 2)), 1));  title("FFT v\_profile without PADDING"); xlim([0, (2 ^ nextpow2(numFrames))]);

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

    fit_options = optimoptions("lsqnonlin", "Display", "off", "Algorithm", "trust-region-reflective");
    try
        [p_fit, ~, residual, exitflag, ~, ~, jacobian] = lsqnonlin(costFunctionHandle, p_init, lb, ub, fit_options);
        [Cn_fit, Dn_fit, center_fit, nu_fit] = modelArgHandle(p_fit);

        ALPHA_N = fitParams.R0 * sqrt(OMEGA_N / nu_fit);
        %          generate_womersley_profile(x, Cn, Dn, center, alpha, width_norm, psf_kernel)
        uWom_fit = generate_womersley_profile(x_coords, Cn_fit, Dn_fit, center_fit, ALPHA_N, width_norm, PSF_KERNEL);

        [Kappa_n, res_mag_RMS, res_phase_RMS, res_phase_RMS_msk, harmonic_SNR_dB] = calculateFitPrecisionMetrics(v_meas, residual, jacobian);

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

        % fitParams.metrics = calculateSymbols(fitParams, RHO_BLOOD, options);
        % Will only replace what is updated and warn for possible new values 
        % not in base
        fitParams.metrics = updateStruct(fitParams.metrics, calculateSymbols(fitParams, RHO_BLOOD, options));

    catch ME 
        warning_s("[WOMERSLEY] Womersley fit failed for %s (%i, %i, %i): %s", name, circleIdx, branchIdx, n_harmonic, ME.message);
        return;
    end

    D_meas_mag_profile = []; 
    
    % Check if d_profile is provided and valid ([Points x 2 x Time])
    if ~(isscalar(d_profile) && isnan(d_profile)) && ~isempty(d_profile) && ndims(d_profile) == 3 && size(d_profile, 3) == numFrames
        
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
    
    % [parabole_fit_systole, parabole_fit_diastole] = analyse_lumen_size(v_profile_good_idx_sav, SYS_IDXS, DIAS_IDXS);
    % 
    % estimated_width.systole = parabole_fit_systole;
    % estimated_width.diastole = parabole_fit_diastole;
    
    % ============================ [ Figures ] ========================== %
    if params.saveFigures && params.json.exportCrossSectionResults.WomersleySectionImage
        hFig = figure("Visible", "off");
        hold on;
    
        title(sprintf('Womersley Fit (%s, %s) for %s (%d, %d) (Harmonic: %d)', options.ModelType, options.NuType, name, circleIdx, branchIdx, n_harmonic), 'Interpreter', 'none');
        
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
    options = optimoptions("lsqnonlin", "Display", "off");
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

function DC_metrics = calculateDCmetrics(geoParams, RHO_BLOOD, FIXED_NU)
    C0 = geoParams.DC_Amp;
    mu_0 = RHO_BLOOD * FIXED_NU; % Viscosité dynamique (Pa.s)
    R0 = geoParams.R0;

    % Calculs analytiques n=0
    K0 = 0.5;
    Q0 = (pi * R0^2 / 2) * C0;
    tau0 = -(2 * mu_0 / R0) * C0;
    G0 = (4 * mu_0 / R0^2) * C0;

    % Stockage des résultats DC
    DC_metrics = struct(...
        "C0", C0, ...
        "K0", K0, ...
        "Q0", Q0, ...
        "tau0", tau0, ...
        "G0", G0, ...
        "R0", R0, ...
        'mu0', mu_0 ...
    );
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
        profile = conv(profile, psf, "same");
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
        Bn = conv(Bn, psf_kernel, "same");
        if Dn ~= 0
            Psin = conv(Psin, psf_kernel, "same");
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
        Bn_profile_final = conv(Bn_profile, psf_kernel, "same");
        Psi_n_profile_final = conv(Psi_n_profile, psf_kernel, "same");
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