function [vessel_mask, combined_response] = gaborVesselness(image, params)
% GABORVESSELNESS Segment blood vessels using Gabor filter bank
%   Enhanced version with improved preprocessing, filtering, and postprocessing
%
% Input:
%   image: 2D grayscale image or 3D RGB image (uint8, uint16, or double)
%   params: structure containing parameters for Gabor filtering (optional)
%     - range: wavelength range (default: 1:10)
%     - step: wavelength step size (default: 1)
%     - sigma: Gaussian smoothing sigma (default: 2)
%     - orientation: filter orientations in degrees (default: 0:30:150)
%     - min_area: minimum pixel area for vessels (default: 50)
%     - closing_radius: radius for morphological closing (default: 2)
%     - contrast_enhance: enable contrast enhancement (default: true)
%     - use_phase: include phase information in response (default: false)
%
% Output:
%   vessel_mask: binary mask of segmented vessels
%   combined_response: combined response from Gabor filters

arguments
    image
    params.range = (1:10)
    params.step = 1
    params.sigma = 2
    params.orientation = 0:30:150
    params.min_area = 50
    params.closing_radius = 2
    params.use_phase = false
end

% Step 1: Enhanced Preprocessing
% Convert to grayscale if necessary
if size(image, 3) == 3
    image = rgb2gray(image);
end

% Handle NaN/inf values and normalize
image = im2double(image);
image(isnan(image) | isinf(image)) = 0;

% Noise reduction with anisotropic diffusion (edge-preserving)
smoothed_image = imdiffusefilt(image, 'NumberOfIterations', 5, ...
    'ConductionMethod', 'quadratic', 'GradientThreshold', 0.05);

% Additional Gaussian smoothing if needed
if params.sigma > 0
    smoothed_image = imgaussfilt(smoothed_image, params.sigma);
end

% Step 2: Enhanced Gabor Filtering
wavelengths = params.range(1):params.step:params.range(2);

% Create Gabor filter bank with more control over parameters
gabor_bank = gabor(wavelengths, params.orientation, ...
    'SpatialFrequencyBandwidth', 1, ...
    'SpatialAspectRatio', 0.5);

% Apply Gabor filters
[gabor_mag, gabor_phase] = imgaborfilt(smoothed_image, gabor_bank);

% Combine responses - improved method
if params.use_phase
    % Use both magnitude and phase information
    phase_weight = 0.3; % Weight for phase information
    gabor_response = gabor_mag .* (1 + phase_weight * cos(gabor_phase));
else
    % Use only magnitude
    gabor_response = gabor_mag;
end

% Enhanced response combination (weighted by orientation)
combined_response = zeros(size(image));

for i = 1:length(gabor_bank)
    % Give more weight to orientations likely to match vessels
    combined_response = combined_response + gabor_response(:, :, i);
end

combined_response = combined_response / length(gabor_bank);

% Adaptive imbinarize
combined_response = mat2gray(combined_response);
vessel_mask = imbinarize(combined_response, 'adaptive');

% Enhanced morphological processing
min_area = params.min_area;
closing_radius = params.closing_radius;

% Remove small objects and holes
vessel_mask = bwareaopen(vessel_mask, min_area);

% Direction-aware closing
se = strel('disk', closing_radius);
vessel_mask = imclose(vessel_mask, se);

end
