function mask = getHolonetprediction(M0, net)
% GETHOLONETPREDICTION Generates a segmentation mask using a U-Net model
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

% Handle network loading
if nargin < 2 || isempty(net)

    if ~isfolder('Models')
        mkdir('Models');
    end

    if ~isfile('Models\uresnet.onnx')
        % Download the model from Hugging Face
        url = 'https://huggingface.co/noTban/uresnet/resolve/main/uresnet.onnx';
        websave('Models\uresnet.onnx', url);
    end

    net = importONNXNetwork('Models\uresnet.onnx');
end

% if nargin<2 || isempty(net)
%     net = importNetworkFromONNX("Models\unet_resnet34.onnx");
% end

[Nx, Ny] = size(M0);

M0 = imresize(rescale(M0), [512, 512]);

input = rescale(M0);

output = predict(net, input);

mask = sigmoid(output) > 0.5;

mask = imresize(mask, [Nx, Ny]);
end
