% ============================== [ SYMBOLS ] ============================ %

function [Kappa_n, res_mag_RMS, res_phase_RMS, res_phase_RMS_msk, harmonic_SNR_dB] = calculateFitPrecisionMetrics(v_meas, residual, jacobian)

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