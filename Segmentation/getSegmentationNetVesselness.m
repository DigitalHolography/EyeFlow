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

model_path = 'Models\' + model_name + '.onnx';

if ~isfile(model_path)
    % Download the model from Hugging Face
    url = 'https://huggingface.co/DigitalHolography/' + model_name + '/resolve/main/' + model_name;
    fprintf(url)
    websave(model_path, url);
end

if canUseGPU
    device = 'gpu';
else
    device = 'cpu';
end

% Import the ONNX network
try
    % Try the newer function first
    net = importNetworkFromONNX(model_path);
catch
    % Fall back to the older function
    warning('off')
    net = importONNXNetwork(model_path);
    warning('on')
end

% if nargin<2 || isempty(net)
%     net = importNetworkFromONNX("Models\unet_resnet34.onnx");
% end

fprintf("    Use %s to segment retinal vessels\n", model_name)

[Nx, Ny] = size(M0);

M0 = imresize(rescale(M0), [512, 512]);

input = rescale(M0);

% First, determine what type of network you have
if isa(net, 'dlnetwork')
    % For dlnetwork objects
    input_dl = dlarray(input, 'SSCB'); % Convert to dlarray
    output_dl = predict(net, input_dl);
    output = extractdata(output_dl);
else
    % For DAGNetwork objects
    output = predict(net, input, 'ExecutionEnvironment', device);
end

mask = sigmoid(output) > 0.5;

mask = imresize(mask, [Nx, Ny], "nearest");
end
