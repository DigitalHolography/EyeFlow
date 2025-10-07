function combinedCrossSectionAnalysis(Q_results_A, Q_results_V, M0_ff_video, sysIdxList)

ToolBox = getGlobalToolBox;

[~, ~, numFrames] = size(M0_ff_video);

fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
t = linspace(0, numFrames / fs, numFrames);
dt = t(2) - t(1);

if ~isempty(sysIdxList)
    sIdx = sysIdxList(1);
    eIdx = sysIdxList(end);
else
    sIdx = 1;
    eIdx = numFrames;
end

Q_cell_A = Q_results_A.Q_cell;
Q_cell_V = Q_results_V.Q_cell;
[numCircles, numBranches_A] = size(Q_cell_A);
[~, numBranches_V] = size(Q_cell_V);

Q_mat_A = nan(numCircles, numBranches_A, numFrames);
Q_mat_V = nan(numCircles, numBranches_V, numFrames);

for cIdx = 1:numCircles

    for bIdx = 1:numBranches_A

        if ~isempty(Q_cell_A{cIdx, bIdx})
            Q_mat_A(cIdx, bIdx, :) = Q_cell_A{cIdx, bIdx};
        end

    end

    for bIdx = 1:numBranches_V

        if ~isempty(Q_cell_V{cIdx, bIdx})
            Q_mat_V(cIdx, bIdx, :) = Q_cell_V{cIdx, bIdx};
        end

    end

end

Q_A = squeeze(mean(Q_results_A.radius_Q, 1)); % Mean over circles
Q_V = squeeze(mean(Q_results_V.radius_Q, 1)); % Mean over circles

numFramesBis = eIdx - sIdx + 1;
tBis = linspace(sIdx / fs, eIdx / fs, numFramesBis);
numCycles = length(sysIdxList) - 1;
cycleSize = numFramesBis / numCycles;

% Volume_A = sum(Q_A) * dt / 60 * 1e-9; % in m^3
% Volume_V = sum(Q_V) * dt / 60 * 1e-9; % in m^3

% Q_ratio = Q_V ./ Q_A;
Q_diff = Q_V - Q_A;
Q_diff = Q_diff - mean(Q_diff(sIdx:eIdx));

% Detrend signals (remove DC offset)
A_detrended = Q_A - mean(Q_A);
V_detrended = Q_V - mean(Q_V);

% Normalize signals (optional but improves cross-correlation robustness)
A_detrended = A_detrended / std(A_detrended);
V_detrended = V_detrended / std(V_detrended);

% Cross-correlation with normalization
[corr_vals, lags] = xcorr(A_detrended, -V_detrended, 'coeff');
[max_corr, max_idx] = max(corr_vals);
time_lag = lags(max_idx) / fs; % Convert lag index to seconds

%% Transfer Function analysis

% Compute FFTs
Q_A_FT = fft(Q_A);
Q_V_FT = fft(Q_V);

F_TRANS = Q_V_FT ./ Q_A_FT;
freqs = linspace(-fs / 2, fs / 2, numel(F_TRANS));
figure("Visible", "off", "Color", 'w');
semilogy(freqs, fftshift(abs(F_TRANS)), '-k', 'LineWidth', 2);
axis tight;
xlabel('Freq (Hz)'); ylabel('transfer function');
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
grid on;
box on;
set(gca, 'LineWidth', 2);
xlim([0 10])

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_Transfer_function_BVR_AV_mod.png", ToolBox.folder_name)))

figure("Visible", "off", "Color", 'w');
plot(freqs, fftshift(angle(F_TRANS)), '-k', 'LineWidth', 2);
axis tight;
xlabel('Freq (Hz)'); ylabel('transfer function angle');
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
grid on;
box on;
set(gca, 'LineWidth', 2);
xlim([0 10])

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_Transfer_function_BVR_AV_phase.png", ToolBox.folder_name)))

ToolBox.Output.Signals.add('TransFunctionModLog10', fftshift(abs(log10(F_TRANS))), 'log10', freqs, 'Hz');
ToolBox.Output.Signals.add('TransFunctionPhaseDegrees', fftshift(180 / pi * angle((F_TRANS))), 'deg', freqs, 'Hz');

instant_dV = detrend(cumsum(Q_diff(sIdx:eIdx))) / 60 * dt;
[peaks, peaks_idx] = findpeaks(instant_dV, 'MinPeakDistance', cycleSize * 0.8);
[troughs, troughs_idx] = findpeaks(-instant_dV, 'MinPeakDistance', cycleSize * 0.8);
peaks = [instant_dV(1) peaks instant_dV(end)];
peaks_idx = [1 peaks_idx numFramesBis];
sys_mean = mean(peaks);
dias_mean = -mean(troughs);

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

%% Figures

% Figure 0 - Flow Rates
figure("Visible", "off", "Color", 'w');
hold on
plot(t, Q_A, 'r', 'LineWidth', 2);
% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, axP(4)]);
xlabel('Time (s)'); ylabel('Flow Rate (µL/min)');
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])

grid on;
box on;
set(gca, 'LineWidth', 2);
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_flowRate_Artery.png", ToolBox.folder_name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_flowRate_Artery.eps", ToolBox.folder_name)))

figure("Visible", "off", "Color", 'w');
hold on
plot(t, Q_V, 'b', 'LineWidth', 2);
% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, axP(4)]);
xlabel('Time (s)'); ylabel('Flow Rate (µL/min)');
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])

grid on;
box on;
set(gca, 'LineWidth', 2);
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_flowRate_Vein.png", ToolBox.folder_name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_flowRate_Vein.eps", ToolBox.folder_name)))

% Figure 1 - Flow Ratio
figure("Visible", "off", "Color", 'w');
hold on;
plot(t, Q_diff, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Flow Rate (µL/min)');

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);

grid on;
box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_flowRate_Diff.png", ToolBox.folder_name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_flowRate_Diff.eps", ToolBox.folder_name)))

% Figure 2 - Signals and cross-correlation
figure("Visible", "off", "Color", 'w');
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
plot(lags / fs, corr_vals, 'k', 'LineWidth', 1.5);
hold on;
plot(time_lag, max_corr, 'ro', 'MarkerSize', 10);
grid on;
axis tight;
xlabel('Lag (s)'); ylabel('Cross-Correlation');
title(['Peak Lag: ', num2str(time_lag, '%.3f'), ' s | Corr: ', num2str(max_corr, '%.2f')]);
box on;
set(gca, 'LineWidth', 2);

% Figure 4 - Cumsum Flow Ratio
figure("Visible", "off", "Color", 'w');
hold on;
plot(tBis, instant_dV, 'k-', 'LineWidth', 2);
scatter(tBis(peaks_idx), peaks, 'red', 'filled', 'o')
scatter(tBis(troughs_idx), -troughs, 'blue', 'filled', 'o')
yline(sys_mean, '--', sprintf("%.2f µL", sys_mean), 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
yline(dias_mean, '--', sprintf("%.2f µL", dias_mean), 'Color', [0.5 0.5 0.5], 'LineWidth', 2)
xlabel('Time (s)');
ylabel('Volume (µL)');
grid on;

% Format plot
axis padded;
axP = axis;
axis([0, numFrames / fs, axP(3), 1.07 * axP(4)]);

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_flowRate_dV.png", ToolBox.folder_name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_flowRate_dV.eps", ToolBox.folder_name)))

close all

end
