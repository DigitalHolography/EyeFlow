function [reconstructed_signal, harmonics, fundamental_freq] = extract_and_reconstruct_harmonics(signal, num_harmonics, fs)
% EXTRACT_AND_RECONSTRUCT_HARMONICS Extract harmonics and reconstruct signal
%
% Inputs:
%   signal - input signal with periodic pattern
%   num_harmonics - number of harmonics to use for reconstruction
%   fs - sampling frequency (Hz)
%
% Outputs:
%   reconstructed_signal - signal reconstructed from harmonics
%   harmonics - structure containing harmonic information
%   fundamental_freq - estimated fundamental frequency (Hz)

    % Ensure signal is a column vector
    signal = signal(:);
    N = length(signal);
    t = (0:N-1)' / fs;
    
    % Estimate fundamental frequency using autocorrelation
    fundamental_freq = estimate_fundamental_frequency(signal, fs);
    
    % Perform FFT
    Y = fft(signal);
    f = fft_freq_vector(fs, N)';
    
    % Find harmonic frequencies and their magnitudes/phases
    harmonics = find_harmonics(Y, f, fundamental_freq, num_harmonics, fs);
    
    % Reconstruct signal from harmonics
    reconstructed_signal = reconstruct_from_harmonics(harmonics, t, N);
    
    % Plot results
    plot_results(signal, reconstructed_signal, harmonics, t, fs);
end

function fundamental_freq = estimate_fundamental_frequency(signal, fs)
% Estimate fundamental frequency using autocorrelation
    
    % Remove DC component
    signal = signal - mean(signal);
    
    % Autocorrelation
    [autocorr, lags] = xcorr(signal, 'coeff');
    autocorr = autocorr(lags >= 0);
    lags = lags(lags >= 0);
    
    % Find peaks in autocorrelation (excluding zero lag)
    [peaks, locs] = findpeaks(autocorr(2:end));
    locs = locs + 1; % Adjust for skipping zero lag
    
    if isempty(locs)
        fundamental_freq = fs / length(signal);
        warning('No clear periodicity found. Using default frequency.');
        return;
    end
    
    % Fundamental period (in samples)
    fundamental_period = lags(locs(1));
    fundamental_freq = fs / fundamental_period;
end

function harmonics = find_harmonics(Y, f, fundamental_freq, num_harmonics, fs)
% Extract harmonic information from FFT
    
    N = length(Y);
    harmonics.frequencies = zeros(num_harmonics, 1);
    harmonics.magnitudes = zeros(num_harmonics, 1);
    harmonics.phases = zeros(num_harmonics, 1);
    
    for k = 1:num_harmonics
        harmonic_freq = k * fundamental_freq;
        
        % Find the closest frequency bin
        [~, idx] = min(abs(f - harmonic_freq));
        
        % Ensure we don't exceed Nyquist frequency
        if harmonic_freq > fs/2
            harmonics.frequencies(k) = NaN;
            harmonics.magnitudes(k) = 0;
            harmonics.phases(k) = 0;
            continue;
        end
        
        harmonics.frequencies(k) = f(idx);
        harmonics.magnitudes(k) = 2 * abs(Y(idx)) / N; % Scale for single-sided spectrum
        harmonics.phases(k) = angle(Y(idx));
    end
    
    % Remove harmonics beyond Nyquist frequency
    valid_idx = ~isnan(harmonics.frequencies);
    harmonics.frequencies = harmonics.frequencies(valid_idx);
    harmonics.magnitudes = harmonics.magnitudes(valid_idx);
    harmonics.phases = harmonics.phases(valid_idx);
end

function reconstructed_signal = reconstruct_from_harmonics(harmonics, t, N)
% Reconstruct signal from harmonic components
    
    reconstructed_signal = zeros(N, 1);
    
    for k = 1:length(harmonics.frequencies)
        if harmonics.magnitudes(k) > 0
            reconstructed_signal = reconstructed_signal + ...
                harmonics.magnitudes(k) * cos(2 * pi * harmonics.frequencies(k) * t + harmonics.phases(k));
        end
    end
end

function plot_results(original_signal, reconstructed_signal, harmonics, t, fs)
% Plot original vs reconstructed signal and harmonic spectrum
    
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot 1: Original vs Reconstructed Signal
    subplot(2, 2, [1, 2]);
    plot(t, original_signal, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original');
    hold on;
    plot(t, reconstructed_signal, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reconstructed');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Original vs Reconstructed Signal');
    legend('show');
    grid on;
    
    % Plot 2: Harmonic Magnitude Spectrum
    subplot(2, 2, 3);
    stem(harmonics.frequencies, harmonics.magnitudes, 'filled', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Harmonic Magnitude Spectrum');
    grid on;
    
    % Plot 3: Error
    subplot(2, 2, 4);
    error_signal = original_signal - reconstructed_signal;
    plot(t, error_signal, 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Error');
    title('Reconstruction Error');
    grid on;
    
    % Display information
    fprintf('Fundamental frequency: %.2f Hz\n', harmonics.frequencies(1));
    fprintf('Number of harmonics used: %d\n', length(harmonics.frequencies));
    fprintf('Reconstruction RMSE: %.6f\n', rms(error_signal));
end