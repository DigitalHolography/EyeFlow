function [signalPlotMean, signalPlotFrames] = signalPlot(signal, stdsignal, opt)

arguments
    signal
    stdsignal
    opt.Visible = "off"
end

% Get global toolbox settings
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

T = ToolBox.stride / ToolBox.fs / 1000;
numFrames = length(signal);
fullTime = linspace(0, numFrames * T, numFrames);

Color_std = [0.7, 0.7, 0.7];
mean_signal = mean(signal);
yAx = [min(-1, min(signal)), max(signal) * 1.07];
axss = [0, numFrames * T, yAx(1), yAx(2)];

% Initialize signal plot with optimized properties
signalPlot = figure('Name', 'Signal Plot', 'Color', 'w', ...
    'Visible', opt.Visible, ...
    'Position', [200 200 600 300]);

% Precompute plot elements
curve1 = signal + stdsignal;
curve2 = signal - stdsignal;
tmp_fullTime = [fullTime, fliplr(fullTime)];
inBetween = [curve1, fliplr(curve2)];

% Create plot with all elements at once (faster than incremental additions)
fill(tmp_fullTime, inBetween, Color_std, 'EdgeColor', 'none');
hold on;
plot(fullTime, curve1, 'Color', Color_std, 'LineWidth', 2);
plot(fullTime, curve2, 'Color', Color_std, 'LineWidth', 2);
plot(fullTime, signal, '-k', 'LineWidth', 2);
yline(mean_signal, '--k', 'LineWidth', 2);

% Set axis properties
axis(axss);

% Formatting
ylabel('Volume Rate (µL/min)');
xlabel('Time (s)');
title(sprintf("Average Blood Volume Rate : %.0f %s", round(mean_signal), 'µL/min'));

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

fontsize(gca, 14, 'points')

signalPlotMean = frame2im(getframe(signalPlot));

% Preallocate signal plot data array
initialFrame = frame2im(getframe(signalPlot));
signalPlotFrames = zeros([size(initialFrame), numFrames], 'uint8');

if exportVideos

    % Create mask coordinates
    maskY = [yAx(1), yAx(2), yAx(2), yAx(1)];

    % Store plot data needed for recreation
    plotData = struct( ...
        'fullTime', fullTime, ...
        'curve1', curve1, ...
        'curve2', curve2, ...
        'signal', signal, ...
        'mean_signal', mean_signal, ...
        'Color_std', Color_std, ...
        'axss', axss, ...
        'T', T, ...
        'numFrames', numFrames, ...
        'maskY', maskY);

    % Parallel frame processing
    parfor frameIdx = 1:numFrames
        % Create temporary figure
        tempFig = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [200 200 600 300]);

        % Recreate the plot directly (faster than copyobj)
        fill([plotData.fullTime, fliplr(plotData.fullTime)], ...
            [plotData.curve1, fliplr(plotData.curve2)], ...
            plotData.Color_std, 'EdgeColor', 'none');
        hold on;
        plot(plotData.fullTime, plotData.curve1, 'Color', plotData.Color_std, 'LineWidth', 2);
        plot(plotData.fullTime, plotData.curve2, 'Color', plotData.Color_std, 'LineWidth', 2);
        plot(plotData.fullTime, plotData.signal, '-k', 'LineWidth', 2);
        yline(plotData.mean_signal, '--k', 'LineWidth', 2);

        % Add frame-specific mask
        patch([frameIdx * plotData.T, frameIdx * plotData.T, ...
                   plotData.numFrames * plotData.T, plotData.numFrames * plotData.T], ...
            plotData.maskY, 'w', 'EdgeColor', 'none', 'FaceAlpha', 1);

        % Set axis properties
        axis(plotData.axss);
        ylabel('Volume Rate (µL/min)');
        xlabel('Time (s)');
        title(sprintf("Average Blood Volume Rate : %.0f %s", round(plotData.mean_signal), 'µL/min'));
        box on;
        set(gca, 'Linewidth', 2, 'Layer', 'top');
        pbaspect([2.5, 1, 1]);
        fontsize(gca, 14, 'points')

        % Capture frame
        signalPlotFrames(:, :, :, frameIdx) = frame2im(getframe(tempFig));
        close(tempFig);
    end

end

end
