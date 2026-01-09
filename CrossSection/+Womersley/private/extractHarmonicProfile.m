% ======================================================================= %
    
function v_meas = extractHarmonicProfile(v_profile_ft, f_vector, base_frequency, n_harmonic, dimToMean)
    % extractHarmonicProfile Extracts the complex velocity profile for a specific harmonic.
    %
    %   v_profile_ft   - The full, FFT-shifted Fourier spectrum of the velocity data.
    %   f_vector       - The corresponding frequency vector for the spectrum.
    %   base_frequency - The fundamental cardiac frequency (e.g., heart rate in Hz).
    %   n_harmonic     - The integer of the harmonic to extract (e.g., 1, 2, 3...).
    %
    % OUTPUT:
    %   v_meas         - The complex velocity profile (a column vector) averaged
    %                    over the frequency band of the specified harmonic.

    arguments
        v_profile_ft, f_vector, base_frequency, n_harmonic
        dimToMean = 2
    end

    target_frequency = n_harmonic * base_frequency;

    freq_resolution = f_vector(2) - f_vector(1);
    margin_hz = freq_resolution * 1.5;
    
    harmonic_indices = find(abs(f_vector - target_frequency) <= margin_hz);

    if isempty(harmonic_indices)
        warning("[WOMERSLEY] Harmonic %d (%.2f Hz) not found in frequency band. Using closest single peak instead.", n_harmonic, target_frequency);
        [~, closest_idx] = min(abs(f_vector - target_frequency));
        harmonic_indices = closest_idx;
    end

    v_meas = mean(v_profile_ft(:, harmonic_indices), dimToMean);
end