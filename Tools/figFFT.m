function figFFT(x)

% Compute the FFT
X = fft(x);

% Compute the magnitude of the FFT (absolute value)
X_mag = log(abs(X));

% Create frequency axis
N = length(x);
fs = 33*10^3/512;  % Sampling frequency (normalized to 1 if unknown)
f = (0:N-1)*(fs/N);  % Frequency vector

% Plot the FFT magnitude
figure;
plot(f, X_mag, 'k', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT of the Signal');
pbaspect([1.618 1 1]);
box on;
set(gca, 'LineWidth', 2)
fontsize(gca, 14, "points") ;

% Adjust axes and labels
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);
hold off;

% X_shifted = log(fftshift(X));
% f_shifted = (-N/2:N/2-1)*(fs/N);  % Centered frequency vector
% figure;
% plot(f_shifted, abs(X_shifted), 'k', 'LineWidth', 2);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Centered FFT of the Signal');
% pbaspect([1.618 1 1]);
% box on;
% set(gca, 'LineWidth', 2)
% fontsize(gca, 14, "points") ;
% 
% % Adjust axes and labels
% axis padded;
% axP = axis;
% axis tight;
% axT = axis;
% axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);
% hold off;
end