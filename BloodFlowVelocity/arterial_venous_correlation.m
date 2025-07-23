function [time_lag, max_corr, lags, corr_vals] = arterial_venous_correlation(v_artery_signal, v_vein_signal)
%ARTERIAL_VENOUS_CORRELATION Computes cross-correlation biomarkers between arterial and venous signals.
%   [time_lag, max_corr, lags, corr_vals] = arterial_venous_correlation(v_artery_signal, v_vein_signal, fs, t)
%
%   Inputs:
%       v_artery_signal : Arterial signal (row/column vector)
%       v_vein_signal   : Venous signal (row/column vector, same length as artery)
%
%   Outputs:
%       time_lag        : Time delay (s) at max cross-correlation
%       max_corr        : Peak cross-correlation value (normalized)
%       lags            : Time lags (s) for cross-correlation
%       corr_vals       : Cross-correlation values

ToolBox = getGlobalToolBox();
fs = ToolBox.fs;
numFrames = length(v_artery_signal);
t = linspace(0, numFrames * ToolBox.stride / fs / 1000, numFrames);

% Input validation
if nargin < 2
    error('Insufficient inputs. Requires at least artery and vein signals.');
end

if length(v_artery_signal) ~= length(v_vein_signal)
    error('Artery and vein signals must be of equal length.');
end

if ~exist('t', 'var') || isempty(t)
    t = (0:length(v_artery_signal) - 1) / fs; % Default time vector
end

% Detrend signals (remove DC offset)
A = v_artery_signal - mean(v_artery_signal);
V = v_vein_signal - mean(v_vein_signal);

% Normalize signals (optional but improves cross-correlation robustness)
A = A / std(A);
V = V / std(V);

% Cross-correlation with normalization
[corr_vals, lags] = xcorr(A, V, 'coeff');
[max_corr, max_idx] = max(corr_vals);
time_lag = lags(max_idx) / fs; % Convert lag index to seconds

% Plot results
figure("Visible", "off");
subplot(2, 1, 1);
hold on
plot(t, A, 'r', 'LineWidth', 2);
plot(t, -V, 'b', 'LineWidth', 2);
axis tight;
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Arterial vs. Venous Signals (Detrended)');
box on;
set(gca, 'LineWidth', 2);

subplot(2, 1, 2);
plot(lags / fs, corr_vals, 'k', 'LineWidth', 1.5);
hold on;
plot(time_lag, max_corr, 'ro', 'MarkerSize', 10);
axis tight;
grid on;
xlabel('Lag (s)'); ylabel('Cross-Correlation');
title(['Peak Lag: ', num2str(time_lag, '%.3f'), ' s | Corr: ', num2str(max_corr, '%.2f')]);
box on;
set(gca, 'LineWidth', 2);

exportgraphics(gcf, fullfile(ToolBox.path_png, ...
    sprintf("%s_arterial_venous_correlation.png", ToolBox.folder_name)))

ToolBox.Outputs.add('PhaseDelay', time_lag, 's', NaN);

close all

end
