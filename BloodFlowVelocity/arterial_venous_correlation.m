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
params = ToolBox.getParams;
saveFigures = params.saveFigures;
fs = ToolBox.fs * 1000 / ToolBox.stride;
numFrames = length(v_artery_signal);
t = ToolBox.Cache.t;

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
[~, max_idx] = max(cross_corr_AV);
max_corr = cross_corr_AV(max_idx);
time_lag = lags(max_idx) / fs; % Convert lag index to seconds
lags_t = lags / fs; % Convert lags to seconds

% MSC calculation
f0 = ToolBox.Cache.HeartBeatFFT;
win_length = 64; % Choose appropriate length for your data
[MSC, f] = mscohere(A, V, hamming(win_length), [], [], fs);
df = 0.3;
valid_indx = (f < (f0 + df)) & (f > (f0 - df));
Gamma_0 = sum(MSC(valid_indx)) ./ sum(valid_indx);

if saveFigures
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

end

ToolBox.Output.DimOut.add('ArteryVeinPhaseDelay', time_lag, ["example", "desc"], 's', NaN);

close all

% Transfer Function analysis

% Compute FFTs
v_a_FT = fft(v_artery_signal, 10 * numFrames);
v_v_FT = fft(v_vein_signal, 10 * numFrames);

F_TRANS = v_v_FT ./ v_a_FT;
freqs = linspace(-fs / 2, fs / 2, numel(F_TRANS));

if saveFigures
    % Plot Transfer Function Magnitude
    figure("Visible", "off", "Color", 'w');
    semilogy(freqs, fftshift(abs(F_TRANS)), '-k', 'LineWidth', 2);
    axis tight;
    xlabel('Freq (Hz)'); ylabel('transfer function');
    set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
    grid on;
    box on;
    set(gca, 'LineWidth', 2);
    xlim([0 10])

    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_Transfer_function_Velocity_AV_mod.png", ToolBox.folder_name)))

    % Plot Phase
    figure("Visible", "off", "Color", 'w');
    plot(freqs, fftshift(angle(F_TRANS)), '-k', 'LineWidth', 2);
    axis tight;
    xlabel('Freq (Hz)'); ylabel('transfer function angle');
    set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
    grid on;
    box on;
    set(gca, 'LineWidth', 2);
    xlim([0 10])

    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_Transfer_function_Velocity_AV_phase.png", ToolBox.folder_name)))
end

ToolBox.Output.Signals.add('TransFunctionModLog10', fftshift(abs(log10(F_TRANS))), 'log10', freqs, 'Hz');
ToolBox.Output.Signals.add('TransFunctionPhaseDegrees', fftshift(180 / pi * angle((F_TRANS))), 'deg', freqs, 'Hz');

% figure;
% subplot(2,1,1); plot(f, 20*log10(abs(Z(1:nfft/2+1))));
% title('Transfer Function Magnitude'); xlabel('Frequency (Hz)'); ylabel('dB');
% subplot(2,1,2); plot(f, angle(Z(1:nfft/2+1)));
% title('Phase Response'); xlabel('Frequency (Hz)'); ylabel('Radians');
%
% Pxx = s_A .* conj(s_A) / nfft;
% Pyy = s_V .* conj(s_V) / nfft;
% Pxy = s_V .* conj(s_A) / nfft;
% Cxy = abs(Pxy).^2 ./ (Pxx .* Pyy);

end
