function mask = getSegmentationNetVesselness(M0, net)
% getSegmentationNetVesselness Generates a segmentation mask using a U-Net model
%   Inputs:
%       M0 - Input image (2D matrix)
%       net - Optional pre-loaded neural network
%   Output:
%       mask - Binary segmentation mask

% Validate input
if nargin < 1 || isempty(M0)
    error('Input image M0 is required');
end

if ~ismatrix(M0)
    error('Input M0 must be a 2D matrix');
end

model_name = "iternet5_vesselness";

if ~isfolder('Models')
    mkdir('Models');
end

if ~isfile('Models\' + model_name +'.onnx')
    % Download the model from Hugging Face
    url = 'https://huggingface.co/DigitalHolography/' + model_name + '/resolve/main/' + model_name;
    fprintf(url)
    websave('Models\' + model_name +'.onnx', url);
end

% Import the ONNX network
warning('off')
net = importONNXNetwork('Models\' + model_name +'.onnx');
warning('on')

% if nargin<2 || isempty(net)
%     net = importNetworkFromONNX("Models\unet_resnet34.onnx");
% end

fprintf("Use " + model_name + " to segment retinal vessels\n")

[Nx, Ny] = size(M0);

M0 = imresize(rescale(M0), [512, 512]);

input = rescale(M0);

output = predict(net, input);

mask = sigmoid(output) > 0.5;

mask = imresize(mask, [Nx, Ny], "nearest");
end
