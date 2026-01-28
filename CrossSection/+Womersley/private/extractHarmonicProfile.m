% ======================================================================= %
    
function v_meas = extractHarmonicProfile(velocityProfileFFT, harmonicNum)
    % extractHarmonicProfile Extracts the complex velocity profile for a specific harmonic.
    %
    %   velocityProfileFFT   - The full, FFT-shifted Fourier spectrum of the velocity data.
    %   f_vector       - The corresponding frequency vector for the spectrum.
    %   base_frequency - The fundamental cardiac frequency (e.g., heart rate in Hz).
    %   harmonicNum     - The integer of the harmonic to extract (e.g., 1, 2, 3...).
    %
    % OUTPUT:
    %   v_meas         - The complex velocity profile (a column vector) averaged
    %                    over the frequency band of the specified harmonic.

    arguments
        velocityProfileFFT
        harmonicNum
    end

    % This means that the velocityProfileFFT is averaged on its cycles
    % And then interpolated to the closest power of 2 (upper)

    % This implies that the FFT is pure, the index 1 is the mean, 2 is
    % cardiac frequency, and so on.

    v_meas = velocityProfileFFT(:, harmonicNum + 1) * 2;

    %{
    target_frequency = harmonicNum * base_frequency;
    
    freq_resolution = f_vector(2) - f_vector(1);
    margin_hz = freq_resolution * 1.5;
    
    harmonic_indices = find(abs(f_vector - target_frequency) <= margin_hz);
    
    if isempty(harmonic_indices)
        warning("[WOMERSLEY] Harmonic %d (%.2f Hz) not found in frequency band. Using closest single peak instead.", harmonicNum, target_frequency);
        [~, closest_idx] = min(abs(f_vector - target_frequency));
        harmonic_indices = closest_idx;
    end
    v_meas = mean(velocityProfileFFT(:, harmonic_indices), options.dimToMean);
    %}
end