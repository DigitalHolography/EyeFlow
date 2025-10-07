function [tau_RC, R_rel, C_rel] = arterial_venous_delay(v_artery, v_vein)

ToolBox = getGlobalToolBox;

v_artery = double(v_artery);
v_vein = double(v_vein);

v_artery_n = rescale(v_artery);
v_vein_n = rescale(v_vein);

numInterp = length(v_artery);

[~, amin] = min(v_vein); % possibly take max(v_artery) instead if not relyable
vein_shift = amin / numInterp * (1 / (ToolBox.Output.HeartBeat.value / 60)); % in seconds

t = (1:numInterp);

tau = fit_tau(t, v_artery_n, v_vein_n, 50, amin);

% ODE definition: dvvein/dt = (v_artery(t) - v_vein) / tau
vvein = vein_solution_conv(t, v_artery_n, 40, 0);

ti = linspace(0, 1 / (ToolBox.Output.HeartBeat.value / 60), numInterp);

% Create figure
hFig = figure('Visible', 'on', 'Color', 'w');
plot(ti, v_artery_n, 'r-', 'Linewidth', 2), hold on
plot(ti, circshift(vvein, amin), 'b--', 'Linewidth', 2)
plot(ti, v_vein_n * max(vvein), 'b-', 'Linewidth', 2)
axis tight;

xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Arterio-venous decay fit', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Arterio-venous decay fit', 'FontSize', 14, 'FontWeight', 'bold');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Compute tau_RC in ms
tau_RC = tau / numInterp * (1 / (ToolBox.Output.HeartBeat.value / 60)); % in seconds
tau_ms = tau_RC * 1000; % convert to ms

% Add legend with tau value
legend({'Artery (normalized)', ...
            sprintf('Vein shifted (%.2f ms)', vein_shift * 1000), ...
            sprintf('Vein model fit (\\tau_{RC} = %.2f ms)', tau_ms)});

R_rel = tau / numInterp;
C_rel = numInterp / tau;
tau_delay_ms = amin / numInterp * 1000 * (1 / (ToolBox.Output.HeartBeat.value / 60));
% Add legend with tau value
legend({'Artery (normalized)', ...
            sprintf('Vein model fit (\\tau_{RC} = %.2f ms)', tau_ms), ...
            sprintf('Vein shifted (%.2f ms)', tau_delay_ms)}, ...
    'Location', 'best');

% Save Results
exportgraphics(hFig, fullfile(ToolBox.path_png, ...
    sprintf("%s_ArterialVenous_Delay.png", ToolBox.folder_name)), ...
    'Resolution', 300);

end

function tau_opt = fit_tau(t, v_artery, v_vein_true, tau0, amin)
%FIT_TAU Find optimal tau that minimizes RSE between model and true vein signal
%
%   tau_opt = fit_tau(t, v_artery, v_vein_true, tau0)
%
% Inputs:
%   t            - time vector (monotonic, uniform sampling assumed)
%   v_artery     - arterial signal (vector)
%   v_vein_true  - measured venous signal (vector, already rescaled if needed)
%   tau0         - initial guess for tau (scalar, optional, default = mean(t))
%
% Output:
%   tau_opt      - optimal tau that minimizes residual squared error

if nargin < 4
    tau0 = mean(t);
end

% Ensure column vectors
t = t(:);
v_artery = v_artery(:);
v_vein_true = v_vein_true(:);

v_vein_sol = @(tau_) vein_solution_conv(t, v_artery, tau_, 0); % force v_vein0 = 0

% Define objective: Residual Squared Error
% force the scale to fit only the shape
target_zone = zeros(length(v_vein_true), 1);
target_zone(amin:end) = 1;
fun = @(tau_) sum(((rescale(v_vein_true) - circshift(rescale(v_vein_sol(tau_)), amin)) .^ 2) .* target_zone);

% Optimization (bounded positive tau)
opts = optimset('TolX', 1e-6);
tau_opt = fminsearch(fun, tau0, opts);
end

% ------------------------------------------------------------------------
function v_vein = vein_solution_conv(t, v_artery, tau, v0)
% Analytic convolution solution of ODE: tau dv/dt + v = v_artery
dt = mean(diff(t));
kernel = exp(-t / tau) / tau;
v_vein = dt * conv(v_artery, kernel, 'full');
v_vein = v_vein(1:length(t)) + v0 * exp(-t / tau);
end
