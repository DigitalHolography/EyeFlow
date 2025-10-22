function mask = getSegmentationNetVesselness(M0, net, ~)
% getSegmentationNetVesselness Generates a segmentation mask using a U-Net model
%   Inputs:
%       M0 - Input image (2D matrix)
%       net - Optional pre-loaded neural network
%   Output:
%       mask - Binary segmentation mask

if canUseGPU
    device = 'gpu';
else
    device = 'cpu';
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
