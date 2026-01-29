function [VelocitySignalPerBeat, VelocitySignalPerBeatFFT, VelocitySignalPerBeatBandLimited] = perBeatSignalAnalysis(signal, sysIdxList, bandLimitedSignalHarmonicCount)
    arguments
        signal
        sysIdxList
        bandLimitedSignalHarmonicCount
    end

    numberOfBeats = numel(sysIdxList) - 1;
    N_fft = 2 ^ nextpow2(max(diff(sysIdxList)));
    VelocitySignalPerBeat = NaN(numberOfBeats, N_fft);
    VelocitySignalPerBeatFFT = NaN(numberOfBeats, N_fft);
    VelocitySignalPerBeatBandLimited = NaN(numberOfBeats, N_fft);
    % perform the fft on each cycle
    for beatIdx = 1:numberOfBeats
        beat = signal(sysIdxList(beatIdx):(sysIdxList(beatIdx + 1)));
        beatInterp = interpft(beat, N_fft + 1);
        beatInterp = beatInterp(1:end-1);
        beatFFT = fft(beatInterp, N_fft);

        VelocitySignalPerBeatFFT(beatIdx,:) = beatFFT;
        VelocitySignalPerBeat(beatIdx,:) = beatInterp;

        bandLimitedSpectrum = beatFFT(1:bandLimitedSignalHarmonicCount) * 2;
        bandLimitedSpectrum(1) = beatFFT(1);

        VelocitySignalPerBeatBandLimited(beatIdx,:) = abs(ifft(bandLimitedSpectrum, N_fft));
    end
end
