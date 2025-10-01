function [best_offsets, scores] = match_triangular(signals, period, rise_ratio)
    %MATCH_TRIANGULAR_MATRIX
    %   Evaluate similarity of multiple signals (N x Nt) to a triangular wave
    %
    %   [best_offsets, scores] = match_triangular_matrix(signals, period, rise_ratio)
    %
    %   signals      : N x Nt matrix (N signals, each of length Nt)
    %   period       : expected period of the triangular pattern
    %   rise_ratio   : fraction of period for rising edge (default = 0.2)
    %
    %   best_offsets : N x 1 vector of offsets for best match
    %   scores       : N x 1 vector of similarity scores
    
    if nargin < 3
        rise_ratio = 0.2;
    end
    
    [N, Nt] = size(signals);
    
    % Generate reference triangular wave (same for all signals)
    ref = triangular_wave(Nt, period, rise_ratio);
    ref = (ref - mean(ref)) ./ (std(ref) + eps); % normalize
    
    % Allocate outputs
    best_offsets = zeros(N,1);
    scores = zeros(N,1);
    
    % Process each signal independently
    for k = 1:N
        sig = signals(k, :);
        sig = (sig - mean(sig)) ./ (std(sig) + eps); % normalize
        
        % Cross-correlation with reference
        [corr_vals, lags] = xcorr(sig, ref, 'normalized');
        
        [scores(k), idx] = max(corr_vals);
        best_offsets(k) = lags(idx);
    end
end


function wave = triangular_wave(len, period, rise_ratio)
    %TRIANGULAR_WAVE Generate asymmetric triangular wave
    wave = zeros(1, len);
    for i = 1:len
        pos = mod(i-1, period);
        if pos < rise_ratio * period
            wave(i) = pos / (rise_ratio * period); % steep rise
        else
            wave(i) = 1 - (pos - rise_ratio * period) / ((1 - rise_ratio) * period);
        end
    end
end
