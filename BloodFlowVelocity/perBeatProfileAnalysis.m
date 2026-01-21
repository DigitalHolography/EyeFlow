function [ProfilePerBeat, ProfilePerBeatFFT, ProfilePerBeatBandLimited] = perBeatProfileAnalysis(profile, sysIdxList, bandLimitedSignalHarmonicCount)
    arguments
        profile double
        sysIdxList double
        bandLimitedSignalHarmonicCount double {mustBePositive}
    end

    numberOfBeats = numel(sysIdxList) - 1;
    N_fft = 2^nextpow2(max(diff(sysIdxList)));
    [~, numPixels] = size(profile);
    ProfilePerBeat            = nan(numberOfBeats, N_fft, numPixels);
    ProfilePerBeatFFT         = nan(numberOfBeats, N_fft, numPixels);
    ProfilePerBeatBandLimited = nan(numberOfBeats, N_fft, numPixels);

    for beatIdx = 1:numberOfBeats
        beatProfile = profile(sysIdxList(beatIdx):sysIdxList(beatIdx+1) - 1, :);

        for p = 1:numPixels
            beat = beatProfile(:, p);
            beatInterp = interpft(beat, N_fft);
            beatFFT = fft(beatInterp, N_fft);

            ProfilePerBeat(beatIdx,:,p)    = beatInterp;
            ProfilePerBeatFFT(beatIdx,:,p) = beatFFT;

            bandLimitedSpectrum = beatFFT(1:bandLimitedSignalHarmonicCount) * 2;
            bandLimitedSpectrum(1) = beatFFT(1);

            ProfilePerBeatBandLimited(beatIdx,:,p) = abs(ifft(bandLimitedSpectrum, N_fft));
        end
    end
end