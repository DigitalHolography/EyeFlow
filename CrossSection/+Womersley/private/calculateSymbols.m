function metrics = calculateSymbols(fitParams, rho_blood, options)
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
    metrics.Kn = calculateKFactor(lambda_n);
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