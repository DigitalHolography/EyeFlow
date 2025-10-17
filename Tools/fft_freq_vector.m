function f = fft_freq_vector(Fs,N,pos)
    if nargin < 3
        pos = 0
    end
    if ~pos
        f = ((0:N-1) - floor(N/2)) * (Fs/N);
    else
        f = (0:floor(N/2)) * (Fs/N);
    end
end