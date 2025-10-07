function [time_lag, max_corr, lags, cross_corr_AV] = arterial_venous_correlation(v_artery_signal, v_vein_signal)
%ARTERIAL_VENOUS_CORRELATION Computes cross-correlation biomarkers between arterial and venous signals.
%   [time_lag, max_corr, lags, cross_corr_AV] = arterial_venous_correlation(v_artery_signal, v_vein_signal, fs, t)
%
%   Inputs:
%       v_artery_signal : Arterial signal (row/column vector)
%       v_vein_signal   : Venous signal (row/column vector, same length as artery)
%
%   Outputs:
%       time_lag        : Time delay (s) at max cross-correlation
%       max_corr        : Peak cross-correlation value (normalized)
%       lags            : Time lags (s) for cross-correlation
%       cross_corr_AV       : Cross-correlation values

ToolBox = getGlobalToolBox();
fs = ToolBox.fs * 1000 / ToolBox.stride;
numFrames = length(v_artery_signal);
t = linspace(0, numFrames / fs, numFrames);

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
[cross_corr_AV, lags] = xcorr(A, V, 'coeff');
[~, max_idx] = max(abs(cross_corr_AV));
max_corr = cross_corr_AV(max_idx);
time_lag = lags(max_idx) / fs; % Convert lag index to seconds
lags_t = lags / fs; % Convert lags to seconds

% MSC calculation
f0 = ToolBox.Cache.list.HeartBeatFFT;
win_length = 64; % Choose appropriate length for your data
[MSC, f] = mscohere(A, V, hamming(win_length), [], [], fs);
df = 0.3;
valid_indx = (f < (f0 + df)) & (f > (f0 - df));
Gamma_0 = sum(MSC(valid_indx)) ./ sum(valid_indx);

figure("Visible", "off", "Color", 'w');
plot(f, MSC, '-k', 'LineWidth', 2)
xline(f0, 'k--', sprintf("%0.2f Hz", f0), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom')
dim = [.6 .5 .3 .3];
str = sprintf("$\\Gamma_0 = %.2f$", Gamma_0);
annotation('textbox', dim, 'String', str, ...
    'FitBoxToText', 'on', 'Interpreter', 'latex', 'FontSize', 14)
axis padded;
xlim([f(1) 10])
xlabel('Frequency (Hz)'); ylabel('Magnitude-squared coherence');
box on;
set(gca, 'LineWidth', 2);

exportgraphics(gcf, fullfile(ToolBox.path_png, ...
    sprintf("%s_arterial_venous_msc.png", ToolBox.folder_name)))

% Plot results
figure("Visible", "off");
subplot(2, 1, 1);
hold on
plot(t, A, 'r', 'LineWidth', 2);
plot(t, -V, 'b', 'LineWidth', 2);
axis tight;
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
box on;
set(gca, 'LineWidth', 2);

subplot(2, 1, 2);
plot(lags_t, cross_corr_AV, 'k', 'LineWidth', 1.5);
hold on;
plot(time_lag, max_corr, 'ro', 'MarkerSize', 10);
axis tight;
grid on;
xlabel('Lag (s)'); ylabel('Cross-Correlation');
legend({sprintf("Peak Lag: %.3f s", time_lag), ...
            sprintf("Peak Corr: %.2f", max_corr)}, 'Location', 'Best');
box on;
set(gca, 'LineWidth', 2);

exportgraphics(gcf, fullfile(ToolBox.path_png, ...
    sprintf("%s_arterial_venous_correlation.png", ToolBox.folder_name)))

% Plot results
figure("Visible", "off");
hold on
plot(t, A, 'r', 'LineWidth', 2);
plot(t, -V, 'b', 'LineWidth', 2);
axis padded;
xlim([0 t(end)])
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
box on;
pbaspect([1.618, 1, 1]);
set(gca, 'LineWidth', 2);

exportgraphics(gcf, fullfile(ToolBox.path_png, ...
    sprintf("%s_detrended_signals.png", ToolBox.folder_name)))

figure("Visible", "off");
plot(lags_t, cross_corr_AV, 'k', 'LineWidth', 1.5);
hold on;
plot(time_lag, max_corr, 'ro', 'MarkerSize', 10);
axis padded;
xlim([lags_t(1) lags_t(end)])
grid on;
xlabel('Lag (s)'); ylabel('Cross-Correlation');
legend({sprintf("Peak Corr: %.2f", max_corr), ...
            sprintf("Peak Lag: %.3f s", time_lag)});
box on;
pbaspect([1.618, 1, 1]);
set(gca, 'LineWidth', 2);

exportgraphics(gcf, fullfile(ToolBox.path_png, ...
    sprintf("%s_lags.png", ToolBox.folder_name)))

ToolBox.Output.add('PhaseDelay', time_lag, 's', NaN);

close all

end
