function derived = calculateSegmentDerivedMetrics(fitParams, qc_model, geoParams, rho_blood)
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
    if isempty(idx1) % || ~qc_model.harmonic_valid(idx1)
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