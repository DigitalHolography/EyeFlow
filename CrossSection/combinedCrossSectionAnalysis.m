function combinedCrossSectionAnalysis(Q_results_A, Q_results_V, M0_ff_video)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

[numX, numY, numFrames] = size(M0_ff_video);

Fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
t = (0:numFrames - 1) / Fs;
dt = t(2) - t(1);

Q_cell_A = Q_results_A.Q_cell;
Q_cell_V = Q_results_V.Q_cell;
[numCircles, numBranches] = size(Q_cell_A);

Q_mat_A = nan(numCircles, numBranches, numFrames);
Q_mat_V = nan(numCircles, numBranches, numFrames);

for cIdx = 1:numCircles
    for bIdx = 1:numBranches
        if ~isempty(Q_cell_A{cIdx, bIdx})
            Q_mat_A(cIdx, bIdx, :) = Q_cell_A{cIdx, bIdx};
        end
        if ~isempty(Q_cell_V{cIdx, bIdx})
            Q_mat_V(cIdx, bIdx, :) = Q_cell_V{cIdx, bIdx};
        end
    end
end

Q_A = mean(squeeze(sum(Q_mat_A, 2, 'omitnan')), 1);
Q_V = mean(squeeze(sum(Q_mat_V, 2, 'omitnan')), 1);

Volume_A = sum(Q_A) * dt / 60 * 1e-9; % in m^3
Volume_V = sum(Q_V) * dt / 60 * 1e-9; % in m^3

Q_ratio = Q_V ./ Q_A;
Q_diff = Q_V - Q_A;

% Detrend signals (remove DC offset)
A_detrended = Q_A - mean(Q_A);
V_detrended = Q_V - mean(Q_V);

% Normalize signals (optional but improves cross-correlation robustness)
A_detrended = A_detrended / std(A_detrended);
V_detrended = V_detrended / std(V_detrended);

% Cross-correlation with normalization
[corr_vals, lags] = xcorr(A_detrended, -V_detrended, 'coeff');
[max_corr, max_idx] = max(corr_vals);
time_lag = lags(max_idx) / Fs; % Convert lag index to seconds

%% Transfer Function analysis

% Compute FFTs
nfft = 2^nextpow2(numFrames); % Zero-pad to next power of 2
f = Fs*(0:(nfft/2))/nfft; % Frequency vector

M = 16;
L = 11;
g = bartlett(M);
Ndft = 1024;

[s_A, f_A, t_A] = spectrogram(Q_A, g, L, Ndft, Fs);
[s_V, f_V, t_V] = spectrogram(Q_V, g, L, Ndft, Fs);

Z = s_V \ s_A; % s_Q_V * Z = s_Q_A

figure;
subplot(2,1,1); plot(f, 20*log10(abs(Z(1:nfft/2+1)))); 
title('Transfer Function Magnitude'); xlabel('Frequency (Hz)'); ylabel('dB');
subplot(2,1,2); plot(f, angle(Z(1:nfft/2+1))); 
title('Phase Response'); xlabel('Frequency (Hz)'); ylabel('Radians');

Pxx = s_A .* conj(s_A) / nfft;
Pyy = s_V .* conj(s_V) / nfft;
Pxy = s_V .* conj(s_A) / nfft;
Cxy = abs(Pxy).^2 ./ (Pxx .* Pyy);

%% Figures
% Figure 1 - Flow Ratio
figure("Visible", "on", "Color", 'w');
hold on;
plot(t, Q_ratio, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Flow Ratio');
title('Q_{vein} / Q_{artery}');
grid on;

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some case

% Figure 1 - Flow Ratio
figure("Visible", "on", "Color", 'w');
hold on;
plot(t, Q_diff, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Flow Rate (µL/min)');
title('Q_{vein} - Q_{artery}');
grid on;

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

% Figure 2 - Signals and cross-correlation
figure("Visible", "on", "Color", 'w');
subplot(2, 1, 1);
hold on
plot(t, V_detrended, 'b', 'LineWidth', 2);
plot(t, A_detrended, 'r', 'LineWidth', 2);
axis tight;
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Arterial vs. Venous Signals (Detrended)');
box on;
set(gca, 'LineWidth', 2);
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

subplot(2, 1, 2);
plot(lags / Fs, corr_vals, 'k', 'LineWidth', 1.5);
hold on;
plot(time_lag, max_corr, 'ro', 'MarkerSize', 10);
axis tight;
grid on;
xlabel('Lag (s)'); ylabel('Cross-Correlation');
title(['Peak Lag: ', num2str(time_lag, '%.3f'), ' s | Corr: ', num2str(max_corr, '%.2f')]);
box on;
set(gca, 'LineWidth', 2);

% Figure 3 - Transfer function analysis
figure("Visible", "on", "Color", 'w');

% Magnitude plot
subplot(3, 1, 1);
semilogx(f, 20*log10(magnitude), 'b', 'LineWidth', 1.5);
hold on;
if ~isempty(valid_freqs)
    semilogx(valid_freqs, 20*log10(valid_mag), 'r.', 'MarkerSize', 10);
end
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(['Transfer Function (Avg gain: ', num2str(avg_gain, '%.2f'), ' or ', num2str(20*log10(avg_gain), '%.1f'), ' dB)']);
xlim([f(2) Fs/2]);
set(gca, 'LineWidth', 2);

% Phase plot
subplot(3, 1, 2);
semilogx(f, phase, 'b', 'LineWidth', 1.5);
hold on;
if ~isempty(valid_freqs)
    semilogx(valid_freqs, valid_phase, 'r.', 'MarkerSize', 10);
end
grid on;
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
title(['Average phase: ', num2str(avg_phase, '%.1f'), ' degrees']);
xlim([f(2) Fs/2]);
set(gca, 'LineWidth', 2);

% Coherence plot
subplot(3, 1, 3);
semilogx(f, coherence, 'b', 'LineWidth', 1.5);
hold on;
plot([f(2) Fs/2], [0.5 0.5], 'r--');
grid on;
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Magnitude-Squared Coherence');
xlim([f(2) Fs/2]);
ylim([0 1]);
set(gca, 'LineWidth', 2);

end
