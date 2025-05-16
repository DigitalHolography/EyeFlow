function graphSignalStd(signalPlot, signal, stdsignal, numFrames, y_label, x_label, fig_title, unit, opt)
% Plots on an existing graph the signal and its std

arguments
    signalPlot
    signal % Signal
    stdsignal % Uncertainty of the Signal
    numFrames
    y_label
    x_label
    fig_title
    unit
    opt.ylimm double = []
    opt.cropIndx double = 0
    opt.fullTime
    opt.xLines = []
    opt.xLineLabels = {}
    opt.ToolBox = []
end

mean_signal = mean(signal);

if opt.cropIndx > 0
    signal = signal(1:opt.cropIndx);
    stdsignal = stdsignal(1:opt.cropIndx);
end

Color_std = [0.7, 0.7, 0.7];
figure(signalPlot);

if ~isempty(opt.ToolBox)
    ToolBox = opt.ToolBox;
    fullTime = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);
else % in a parfor no ToolBox
    fullTime = opt.fullTime;
end

if isempty(opt.ylimm)
    axss = [fullTime(1), fullTime(end), min(signal), max(signal)];
else
    axss = [fullTime(1), fullTime(end), opt.ylimm];
end

if length(signal) ~= numFrames % for a variable length of the signal
    fullTime = fullTime(1:length(signal));
end

curve1 = signal + stdsignal;
curve2 = signal - stdsignal;
tmp_fullTime = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

fill(tmp_fullTime, inBetween, Color_std);
hold on;
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2);
plot(fullTime, signal, '-k', 'LineWidth', 2);
yline(mean_signal, '--k', 'LineWidth', 2)
hold off;

ylabel(y_label)
xlabel(x_label)
title(sprintf("%s : %.0f %s", fig_title, round(mean_signal), unit))

if ~isempty(opt.xLines)

    for n = 1:length(opt.xLines)
        xline(opt.xLines(n), ':', obj.xLineLabels{n}, LineWidth = opt.LineWidth);
    end

end

if ~isempty(opt.ylimm)
    axis(axss);
else
    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), 0, 1.07 * axP(4)])
end

fontsize(gca, 14, "points");
set(gca, 'Linewidth', 2)
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])

end
