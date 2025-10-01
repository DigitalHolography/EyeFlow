function selected_idx = select_regular_peaks(signals_n, threshold, tolerance)
%SELECT_REGULAR_PEAKS Select signals with regular derivative peaks
%
%   INPUTS:
%       signals_n : (N x Nt) matrix, each row is a normalized signal
%       threshold : minimum derivative peak amplitude
%       tolerance : allowed relative variation in peak intervals (e.g. 0.1 = 10%)
%
%   OUTPUT:
%       selected_signals : subset of signals_n that show regular peaks

ToolBox = getGlobalToolBox;
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);

    [N, Nt] = size(signals_n);
    selected_idx = [];

    for i = 1:N
        sig = signals_n(i,:);

        % Apply bandpass filter
        [b, a] = butter(4, 15 / (fs / 2), 'low');
        filtered_signal = filtfilt(b, a, sig);


        d_sig = diff(filtered_signal);  % derivative approximation

        if max(d_sig) > threshold

            % Find positive peaks in derivative
            [pks, locs] = findpeaks(d_sig, 'MinPeakHeight', threshold);
    
            if numel(locs) > 2
                % Compute intervals between successive peaks
                intervals = diff(locs);
    
                % Check regularity: coefficient of variation (std/mean)
                cv = std(intervals) / mean(intervals);
    
                if cv < tolerance
                    selected_idx = [ selected_idx i ];
                end
            end
        end
    end

    % Keep only selected signals
    % selected_signals = signals_n(selected_idx,:);
end
