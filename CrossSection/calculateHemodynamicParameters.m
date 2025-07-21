function [R, tau_w, mu, uncertainties] = calculateHemodynamicParameters(Q_t, dQ_t, deltaP, r, L, index_start, index_end, N)
% CALCULATEHEMODYNAMICPARAMETERS Computes key hemodynamic parameters
%   [R, tau_w, mu, uncertainties] = calculateHemodynamicParameters(Q, deltaP, r, L)
%   calculates vascular resistance (R), wall shear stress (tau_w), and dynamic
%   viscosity (mu) using Hagen-Poiseuille relationships.
%
% Inputs:
%   Q       - Volumetric flow rate (m^3/s) [mean, std]
%   deltaP  - Pressure drop (Pa) [mean, std]
%   r       - Vessel radius (m) [mean, std]
%   L       - Vessel length (m) [mean, std]
%
% Outputs:
%   R       - Vascular resistance (Pa·s/m^3) [mean, std]
%   tau_w   - Wall shear stress (Pa) [mean, std]
%   mu      - Dynamic viscosity (Pa·s) [mean, std]
%   uncertainties - Structure with relative uncertainties

% Constants
mmHg_to_Pa = 133.322; % Conversion factor
muLmin_to_m3s = 1e-9 / 60;

Q_mean = mean(Q_t(index_start:index_end)) * muLmin_to_m3s;
Q_std = sqrt(mean(dQ_t(index_start:index_end) .^ 2)) * muLmin_to_m3s; % RMS of the uncertainty

if numel(deltaP) == 2
    deltaP_mean = deltaP(1); deltaP_std = deltaP(2);
else
    deltaP_mean = deltaP; deltaP_std = 0;
end

if numel(r) == 2
    r_mean = r(1); r_std = r(2);
else
    r_mean = r; r_std = 0;
end

if numel(L) == 2
    L_mean = L(1); L_std = L(2);
else
    L_mean = L; L_std = 0;
end

% Core Calculations
% Resistance (R = ΔP/Q)
R_mean = deltaP_mean * N / Q_mean;
R_std = R_mean * sqrt((deltaP_std/deltaP_mean)^2 + (Q_std/Q_mean)^2);

% Viscosity (μ = πr⁴ΔP/(8LQ))
mu_mean = (pi * r_mean^4 * deltaP_mean) / (8 * L_mean * Q_mean);
mu_mean = mu_mean * N; % Only one vessel
mu_std = mu_mean * sqrt(...
    (4*r_std/r_mean)^2 + ...
    (deltaP_std/deltaP_mean)^2 + ...
    (L_std/L_mean)^2 + ...
    (Q_std/Q_mean)^2);

% Wall shear stress (τ_w = 4μQ/(πr³))
tau_w_mean = (4 * mu_mean * Q_mean) / (pi * r_mean^3);
tau_w_mean = tau_w_mean / N; % Only one vessel
tau_w_std = tau_w_mean * sqrt(...
    (mu_std/mu_mean)^2 + ...
    (Q_std/Q_mean)^2 + ...
    (3*r_std/r_mean)^2);

% Prepare outputs
R = [R_mean, R_std];
tau_w = [tau_w_mean, tau_w_std];
mu = [mu_mean, mu_std];

uncertainties = struct(...
    'R_rel', R_std/R_mean, ...
    'tau_w_rel', tau_w_std/tau_w_mean, ...
    'mu_rel', mu_std/mu_mean);

% Display results if requested

fprintf('\n--- Hemodynamic Parameters ---\n');
fprintf('Flow Rate (Q): %.2f ± %.2f mm³/s (%.1f%%)\n', ...
    Q_mean * 1e9 , Q_std * 1e9, 100*Q_std/Q_mean);
fprintf('Vascular Resistance (R): %.1f ± %.1f mmHg·s/mm³ (%.1f%%)\n', ...
    R_mean/1e9/mmHg_to_Pa, R_std/1e9/mmHg_to_Pa, 100*uncertainties.R_rel);
fprintf('Dynamic Viscosity (μ): %.1f ± %.1f mPa·s (%.1f%%)\n', ...
    mu_mean*1000, mu_std*1000, 100*uncertainties.mu_rel);
fprintf('Wall Shear Stress (τ_w): %.1f ± %.1f Pa (%.1f%%)\n\n', ...
    tau_w_mean, tau_w_std, 100*uncertainties.tau_w_rel);

end