function spectrum_ploting(SH, f_signal, f_signal_bkg, mask_artery, mask_bkg, fs, f1, f2)
%SPECTRUM_PLOTING Summary of this function goes here
%   Detailed explanation goes here

[~, ~, numFreq] = size(SH);

f_signal_idx = f_signal * numFreq / fs;
f_signal_idx_bkg = f_signal_bkg * numFreq / fs;

SH_artery = squeeze(sum(SH .* mask_artery, [1 2])) / nnz(mask_artery);
SH_bkg = squeeze(sum(SH .* mask_bkg, [1 2])) / nnz(mask_bkg);

I_f = fftshift(log10(SH_artery(round(f_signal_idx))));
I_f_bkg = fftshift(log10(SH_bkg(round(f_signal_idx_bkg))));

axis_x = linspace(-fs / 2, fs / 2, numFreq);

p_artery = plot(axis_x, fftshift(log10(SH_artery)), 'red', 'LineWidth', 1.5, 'DisplayName', 'Arteries');
sclingrange = abs(fftshift(axis_x)) > f1;
yrange = [.99 * log10(min(SH_artery(sclingrange))) 1.01 * log10(max(SH_artery(sclingrange)))];
ylim(yrange);
xlim([-fs / 2 fs / 2]);
rectangle('Position', [-f1 yrange(1) 2 * f1 (yrange(2) - yrange(1))], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [-fs / 2 yrange(1) (fs / 2 - f2) (yrange(2) - yrange(1))], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
rectangle('Position', [f2 yrange(1) (fs / 2 - f2) (yrange(2) - yrange(1))], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
om_RMS_line = line([-f_signal f_signal], [I_f I_f]);
om_RMS_line.Color = 'red';
om_RMS_line.LineStyle = '-';
om_RMS_line.Marker = '|';
om_RMS_line.MarkerSize = 12;
om_RMS_line.LineWidth = 1;
om_RMS_line.Tag = 'f RMS';
text(0, I_f, 'f_{RMS}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

xline(f1, '--')
xline(f2, '--')
xline(-f1, '--')
xline(-f2, '--')
xticks([-f2 -f1 0 f1 f2])
xticklabels({num2str(round(-f2, 1)), num2str(round(-f1, 1)), '0', num2str(round(f1, 1)), num2str(round(f2, 1))})
% annotation('textarrow',omegaAVG_left,omegaAVG_right,'String','omega','FontSize',13,'Linewidth',2)
title('Average spectrum')

% plot(fullTime,fullArterialPulse,'-k', fullTime,fullBackgroundSignal,':k', fullTime, fullVenousSignal, '-.k', 'LineWidth',2) ;
% title('arterial pulse waveform and background signal'); % averaged outside of segmented vessels
fontsize(gca, 12, "points");
xlabel('frequency (kHz)', 'FontSize', 14);
ylabel('log10 S', 'FontSize', 14);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 1);
uistack(p_artery, 'top');
uistack(gca, 'top');

hold on
plot(axis_x, fftshift(log10(SH_bkg)), 'black--', 'LineWidth', 1.5, 'DisplayName', 'Background');
om_RMS_line = line([-f_signal_bkg f_signal_bkg], [I_f_bkg I_f_bkg]);
om_RMS_line.Color = 'black';
om_RMS_line.LineStyle = '-';
om_RMS_line.Marker = '|';
om_RMS_line.MarkerSize = 12;
om_RMS_line.LineWidth = 1;
om_RMS_line.Tag = 'f RMS';
text(0, I_f_bkg, 'f_{RMS background}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
legend('', '', '', '', '', 'Arteries', 'Background');

pbaspect([1.618 1 1])
box on
set(gca, 'LineWidth', 2)

end
