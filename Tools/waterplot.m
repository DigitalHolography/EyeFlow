function waterplot(s, f, t)
% Waterfall plot of spectrogram
waterfall(f, t, log(abs(s)' .^ 2))
set(gca, XDir = "reverse", View = [30 50])
xlabel("Frequency (Hz)")
ylabel("Time (s)")
end
