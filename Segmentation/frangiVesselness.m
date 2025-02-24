function [M0_binary_img] = frangiVesselness(M0_ff_img, name, TB, opt)
    % Input: M0_ff_img - Noisy vessel image (numX x numY)
    % Output: M0_denoised_img - Denoised vessel image

    arguments
        M0_ff_img
        name 
        TB 
        opt.BlackWhite = false
    end

    % Step 1: Preprocessing
    % Convert to grayscale if the image is RGB
    if size(M0_ff_img, 3) == 3
        M0_gray_img = rgb2gray(M0_ff_img);
    else
        M0_gray_img = M0_ff_img;
    end

    % Normalize the image to [0, 1]
    M0_norm_img = rescale(M0_gray_img);

%     % Step 2: Denoising
%     PW_params = Parameters_json(ToolBox.PW_path,ToolBox.PW_param_name);
%     % Apply Gaussian smoothing
%     sigma = PW_params.gauss_filt_size_for_barycenter; % Standard deviation for Gaussian filter
%     M0_gaussian_img = imgaussfilt(M0_norm_img, sigma);
% 
%     % Apply Non-Local Means (NLM) denoising
%     patch_size = 5; % Size of patches for NLM
%     search_window = 11; % Size of search window for NLM
%     M0_nlm_img = imnlmfilt(M0_gaussian_img, 'DegreeOfSmoothing', 0.05, ...
%                            'SearchWindowSize', search_window, ...
%                            'ComparisonWindowSize', patch_size);

    PW_params = TB.getParams();

    % Step 3: Vessel Enhancement
    % Apply Frangi Vesselness Filter
    M0_vesselness_img = FrangiFilter2D(M0_norm_img, ...
        "FrangiScaleRange", PW_params.params.Mask.VesselnessSigmaRange, ...
        "FrangiScaleRatio", PW_params.params.Mask.VesselnessSigmaStep, ...
        "BlackWhite", opt.BlackWhite);

    % Thresholding to segment vessels (optional)
    M0_binary_img = imbinarize(M0_vesselness_img, "adaptive");

    saveImage(M0_vesselness_img, TB, sprintf('%s_frangi_img.png', name), isStep = true)
    saveImage(M0_binary_img, TB, sprintf('%s_frangi_mask.png', name), isStep = true)
end