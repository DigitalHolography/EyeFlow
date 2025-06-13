function [filtered_signal] = svd_filter(signal, window_size, rank)
% SVD_FILTER Applies SVD-based filtering to a 1D signal
%
% Inputs:
%   signal      - Input signal (1D vector)
%   window_size - Size of the Hankel matrix window (should be <= length(signal)/2)
%   rank        - Number of singular values to keep
%
% Output:
%   filtered_signal - Denoised version of the input signal

N = length(signal);
L = window_size; % Number of rows in Hankel matrix
K = N - L + 1; % Number of columns in Hankel matrix

% Create Hankel matrix
H = zeros(L, K);

for i = 1:L

    for j = 1:K
        H(i, j) = signal(i + j - 1);
    end

end

% Perform SVD
[U, S, V] = svd(H, 'econ');

% Truncate to selected rank
U_r = U(:, 1:rank);
S_r = S(1:rank, 1:rank);
V_r = V(:, 1:rank);

% Reconstruct matrix
H_filtered = U_r * S_r * V_r';

% Diagonal averaging (Hankelization)
filtered_signal = zeros(N, 1);

for n = 1:N

    if n <= L && n <= K
        indices = 1:n;
    elseif n > L && n <= K
        indices = 1:L;
    elseif n <= L && n > K
        indices = n - K + 1:n;
    else
        indices = n - K + 1:L;
    end

    antidiag_indices = sub2ind(size(H_filtered), ...
        fliplr(indices), ...
        n - indices + 1);
    filtered_signal(n) = mean(H_filtered(antidiag_indices));
end

end
