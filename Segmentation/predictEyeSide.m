function eyeSide = predictEyeSide(model, M0)
% predictEyeSide - Predicts if a M0 image shows a left or right eye
%   Inputs:
%       model     - The pre-loaded YOLO model
%       M0        - Input image (2D matrix)
%   Output:
%       eyeSide   - string: "left" or "right"

M0 = rescale(imresize(M0, [640, 640]));
rgbM0 = cat(3, M0, M0, M0);
dlInput = dlarray(rgbM0, 'SSCB');

dlOutput = predict(model, dlInput);

scores = extractdata(gather(dlOutput));
scores = scores(:);

[maxScore, maxIndex] = max(scores);

if maxIndex == 1
    eyeSide = "left";
    labelStr = "Left";
    borderColor = [0 200 0];
else
    eyeSide = "right";
    labelStr = "Right";
    borderColor = [200 0 0];
end

% Save the image
Toolbox = getGlobalToolBox;

if Toolbox.getParams.saveFigures
    visImg = uint8(rgbM0 * 255);
    borderWidth = 10;

    visImg(1:borderWidth, :, 1) = borderColor(1);
    visImg(1:borderWidth, :, 2) = borderColor(2);
    visImg(1:borderWidth, :, 3) = borderColor(3);
    visImg(end - borderWidth + 1:end, :, 1) = borderColor(1);
    visImg(end - borderWidth + 1:end, :, 2) = borderColor(2);
    visImg(end - borderWidth + 1:end, :, 3) = borderColor(3);
    visImg(:, 1:borderWidth, 1) = borderColor(1);
    visImg(:, 1:borderWidth, 2) = borderColor(2);
    visImg(:, 1:borderWidth, 3) = borderColor(3);
    visImg(:, end - borderWidth + 1:end, 1) = borderColor(1);
    visImg(:, end - borderWidth + 1:end, 2) = borderColor(2);
    visImg(:, end - borderWidth + 1:end, 3) = borderColor(3);

    precisionText = sprintf('%s: %.2f%%', labelStr, maxScore * 100);

    try
        visImg = drawTextBox( ...
            visImg, ...
            [borderWidth + 5, borderWidth + 5], ...
            precisionText, ...
            24, ...
            borderColor, ...
            1, ...
        'white');
    catch
        warning('insertText not available. Saving image with border only.');
    end

    Toolbox = getGlobalToolBox;
    imwrite(visImg, fullfile(Toolbox.path_png, 'eye_side.png'));
end

Toolbox.Output.add('PredictedEyeSide', eyeSide, h5path = '/Meta/EyePredictedSide');

end
