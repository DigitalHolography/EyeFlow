function [peak_freqs, peaks, phase] = syntheticSpectralAnalysis(signal, t, fs, N)
% syntheticSpectralAnalysis - Perform spectral analysis on a one-cycle blood flow velocity signal
%
% Inputs:
%   signal          - The blood flow velocity signal (1D array)
%   t               - Time vector corresponding to the signal (1D array)
%   fs              - Estimated fundamental frequency (Hz)
%   N               - Number of synthetic replications to increase frequency resolution
% Outputs:
%   peak_freqs     - Frequencies of the detected peaks (harmonics)
%   peaks          - Magnitudes of the detected peaks
%   phase          - Phases of the detected peaks

% Increase resolution with synthetic replication of the averaged cycle
synthetic_signal = repmat(signal, [1 N]); % N times

% Only get the positive components
f = fft_freq_vector(fs, length(synthetic_signal));
f = f(end / 2 + 1:end);
fft_synthetic = fft(synthetic_signal);
fft_synthetic = fft_synthetic(1:end / 2);
fft_synthetic(2:end / 2) = fft_synthetic(2:end / 2) .* 2;

% Get the Magnitude and Phase
fft_mag = abs(fft_synthetic);
fft_mag = fft_mag ./ fft_mag(1) .* mean(synthetic_signal); % Normalised
fft_angle = angle(fft_synthetic);

% Find harmonics
[peaks, peak_locs] = findpeaks(fft_mag, 'MinPeakHeight', 0.001);

% Add the mean component
peaks = [fft_mag(1), peaks];
peak_locs = [1, peak_locs];

peak_freqs = f(peak_locs);
phase = fft_angle(peak_locs);

% Display
%{
n = (1:length(peak_locs));
newSignal = sum(peaks(n) .* cos(t' .* 2 .* pi * peak_freqs(n) + phase(n)), 2)';
newSignalSimple = sum(peaks(1:6) .* cos(t' .* 2 .* pi * peak_freqs(1:6) + phase(1:6)), 2)';

figure,
subplot(4, 1, 1), plot(t, signal),
axis padded
axP = axis;
axis([t(1), t(end), 0, axP(4)])
xlabel("Time (s)"), ylabel("Blood Velocity (mm/s)")

subplot(4, 1, 2), plot(f, fft_mag),
hold on, scatter(peak_freqs, peaks), xlim([0 10])
xlabel("Frequency (Hz)"), ylabel("Magnitude |Y|")

subplot(4, 1, 3), plot(f, fft_angle / pi),
hold on, scatter(peak_freqs, phase / pi), xlim([0 10]), ylim([-1 1])
xlabel("Frequency (Hz)"), ylabel("Phase / \pi")

subplot(4, 1, 4), plot(t, newSignal),
hold on, plot(t, signal),
plot(t, newSignalSimple, ':')
xlabel("Time (s)"), ylabel("Blood Velocity (mm/s)")
axis([t(1), t(end), 0, axP(4)])
%}

end
