function [maskArtery, maskVein] = assert_arteries(video, maskArtery, maskVein)
% This function makes sure maskArtery and Vein are not inverted by choosing the signal with the largest
% derimative max.
% Ensure maskArtery and maskVein are not inverted by comparing the max derivative of their signals

ToolBox = getGlobalToolBox();
fs = ToolBox.fs * 1000 / ToolBox.stride; % Convert to Hz

signal = sum(video .* maskArtery, [1 2], 'omitnan');
signal = signal ./ nnz(maskArtery);

outlier_frames_mask = isoutlier(signal, "movmedian", 5, "ThresholdFactor", 2);
video = interpolateOutlierFrames(video, outlier_frames_mask);

signal = sum(video .* maskArtery, [1 2], 'omitnan');
signal = signal ./ nnz(maskArtery);

% Extract pulse signals for artery and vein masks
pulseArtery = squeeze(sum(video .* maskArtery, [1 2], 'omitnan') / nnz(maskArtery));
outlier_frames_mask = isoutlier(pulseArtery, "movmedian", 5, "ThresholdFactor", 2);
video_artery = interpolateOutlierFrames(video, outlier_frames_mask);
pulseArtery = squeeze(sum(video_artery .* maskArtery, [1 2], 'omitnan') / nnz(maskArtery));

% pulseVein   = squeeze(sum(video .* maskVein,   [1 2], 'omitnan') / nnz(maskVein));
% outlier_frames_mask = isoutlier(pulseVein, "movmedian", 5, "ThresholdFactor", 2);
% video_vein = interpolateOutlierFrames(video, outlier_frames_mask);
pulseVein = squeeze(sum(video_artery .* maskVein, [1 2], 'omitnan') / nnz(maskVein));

% Filter both signals
[b, a] = butter(4, 15 / (fs / 2), 'low');
pulseArtery = filtfilt(b, a, pulseArtery);
pulseVein = filtfilt(b, a, pulseVein);

% Compute derivatives
diffArtery = gradient(pulseArtery);
diffVein = gradient(pulseVein);

% Compare max derivative
if max(diffVein) > max(diffArtery)
    % Swap masks if vein has higher max derivative
    disp("Inverting artery and vein due to the signals analysis.")
    tmp = maskArtery;
    maskArtery = maskVein;
    maskVein = tmp;
end

end
