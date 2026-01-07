% ================================ [ QC ] =============================== %

function qc = calculateWomersleyQC(fitParams, modelType, varargin)
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
