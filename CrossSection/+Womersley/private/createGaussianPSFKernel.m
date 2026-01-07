function psf_kernel = createGaussianPSFKernel(fwhm_um, num_points, cross_section_px, pixel_size_um)
    % Gaussian PSF kernel

    if fwhm_um <= 0
        % If FWHM is zero or negative, return an empty kernel to signify no convolution.
        psf_kernel = [];
        return;
    end

    total_width_um = cross_section_px * pixel_size_um;
    dx_um = total_width_um / num_points; % um per point

    if dx_um == 0
        warning("[WOMERSLEY] Cannot create PSF kernel: The spatial resolution (dx_um) is zero.");
    end
    
    fwhm_in_points = fwhm_um / dx_um;

    sigma_in_points = fwhm_in_points / (2 * sqrt(2 * log(2)));

    kernel_radius = ceil(3 * sigma_in_points);
    kernel_size = 2 * kernel_radius + 1;

    % 'fspecial' creates a normalized kernel that sums to 1
    psf_kernel = fspecial('gaussian', [1, kernel_size], sigma_in_points);
end
  