function mask = getHolonetprediction(M0, net)

if nargin < 2 || isempty(net)
    net = importONNXNetwork("Models\unet_resnet34.onnx");
end

[Nx, Ny] = size(M0);

M0 = imresize(rescale(M0), [512, 512]);

input = rescale(M0);

output = predict(net, input);

mask = sigmoid(output(:, :, 2)) > 0.5;

mask = imresize(mask, [Nx, Ny]);
end
