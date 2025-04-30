function spectrum_video(SH, f_RMS_video, mask_artery, mask_bkg)

ToolBox = getGlobalToolBox;
[numX, numY, numFreq, numFrames] = size(SH);
f1 = ToolBox.f1;
f2 = ToolBox.f2;
fs = ToolBox.fs;

f_signal = squeeze(sum(f_RMS_video .* mask_artery, [1 2]) / nnz(mask_artery));
f_signal_bkg = squeeze(sum(f_RMS_video .* mask_bkg, [1 2]) / nnz(mask_bkg));

mask_artery = logical(imresize(mask_artery, [numX, numY]));
mask_bkg = logical(imresize(mask_bkg, [numX, numY])) & ~mask_artery;

% make video
spectrum_video = zeros(420, 560, 3, numFrames);

f_signal_idx = f_signal * numFreq / fs;
f_signal_idx_bkg = f_signal_bkg * numFreq / fs;

SH_artery = squeeze(sum(SH .* mask_artery, [1 2])) / nnz(mask_artery);
SH_bkg = squeeze(sum(SH .* mask_bkg, [1 2])) / nnz(mask_bkg);

axis_x = linspace(-fs / 2, fs / 2, numFreq);
sclingrange = abs(fftshift(axis_x)) > f1;
ymin = log10(min(SH_artery(sclingrange, :), [], 'all'));
ymax = log10(max(SH_artery(sclingrange, :), [], 'all'));

parfor frameIdx = 1:numFrames
    Ninterp = 10;
    axis_x = linspace(-fs / 2, fs / 2, Ninterp * numFreq);
    fi = figure("Visible", "on", "Color", 'w');

    SH_t = interp1(linspace(-fs / 2, fs / 2, numFreq), SH_artery(:, frameIdx), axis_x);
    SH_bkg_t = interp1(linspace(-fs / 2, fs / 2, numFreq), SH_bkg(:, frameIdx), axis_x);

    I_f = fftshift(log10(SH_t(round(Ninterp * f_signal_idx(frameIdx)))));
    I_f_bkg = fftshift(log10(SH_bkg_t(round(Ninterp * f_signal_idx_bkg(frameIdx)))));

    rectangle('Position', [-f1 ymin 2 * f1 (ymax - ymin)], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
    hold on

    plot(axis_x, fftshift(log10(SH_t)), 'red', 'LineWidth', 1.5, 'DisplayName', 'Arteries');
    ylim([ymin ymax]);
    xlim([-fs / 2 fs / 2]);
    om_RMS_line = line([-f_signal(frameIdx) f_signal(frameIdx)], [I_f I_f]);
    om_RMS_line.Color = 'red';
    om_RMS_line.LineStyle = '-';
    om_RMS_line.Marker = '|';
    om_RMS_line.MarkerSize = 12;
    om_RMS_line.LineWidth = 1;
    om_RMS_line.Tag = 'f RMS';
    text(0, I_f, 'f_{RMS}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

    plot(axis_x, fftshift(log10(SH_bkg_t)), 'black--', 'LineWidth', 1.5, 'DisplayName', 'Background');
    om_RMS_line = line([-f_signal_bkg(frameIdx) f_signal_bkg(frameIdx)], [I_f_bkg I_f_bkg]);
    om_RMS_line.Color = 'black';
    om_RMS_line.LineStyle = '-';
    om_RMS_line.Marker = '|';
    om_RMS_line.MarkerSize = 12;
    om_RMS_line.LineWidth = 1;
    om_RMS_line.Tag = 'f RMS';
    text(0, I_f_bkg, 'f_{RMS background}', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')

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

    fontsize(gca, 14, "points");
    xlabel('frequency (kHz)', 'FontSize', 14);
    ylabel('log10 S', 'FontSize', 14);
    pbaspect([1.618 1 1]);
    box on
    set(gca, 'LineWidth', 2);

    legend('Arteries', '', 'Background');

    frame = getframe(fi);
    spectrum_video(:, :, :, frameIdx) = frame.cdata;
    close(fi)
end

writeVideoOnDisc(spectrum_video, fullfile(ToolBox.path_avi, strcat(ToolBox.main_foldername, '_spectrum_video')));
writeGifOnDisc(mat2gray(spectrum_video), "spectrum_video")
end
