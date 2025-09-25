function obj = VideoNonRigidRegistering(obj)
tic
params = Parameters_json(obj.directory, obj.param_name);

if ~params.json.Preprocess.NonRigidRegisteringFlag
    return
end

fprintf("    - Video Registering Non-Rigidly started...\n");

ref_img = mean(obj.M0_raw_video, 3);
ref_img = flat_field_correction(ref_img, 35, 0, 'gaussianBlur');

v = obj.M0_raw_video;

outDir = fullfile(obj.directory, 'eyeflow', 'nonRigidReg');

if ~isfolder(outDir)
    mkdir(outDir);
end

save_path = fullfile(obj.directory, 'eyeflow', "nonRigidReg");

nFrames = size(v, 3);

stabilized = zeros(size(ref_img, 1), size(ref_img, 2), nFrames);
field = zeros(size(ref_img, 1), size(ref_img, 2), 2, nFrames);

%smoothVideo = imgaussfilt3(obj.M0_raw_video, [0.1 0.1 2]);
parfor k = 1:nFrames
    %get the frame, stabilize it, save it
    tgt = safeConvertFrame(v(:, :, k));
    tgt = flat_field_correction(tgt, 35, 0, 'gaussianBlur');
    [f, s] = diffeomorphicDemon(tgt, ref_img, tgt);

    field(:, :, :, k) = f;
    stabilized(:, :, k) = s;
end

D.stabilized = stabilized;
D.field = field;

saveStabilizedVideoGif(D.stabilized, fullfile(save_path, "stabilized.gif"));
saveAngleAndNormOfDisplacementField(D.field, fullfile(save_path, "angle.gif"), fullfile(save_path, "norm.gif"));

obj.displacementField = D;
fprintf("    - Video Registering Non-Rigidly took: %ds\n", round(toc));
end

% helper functions

function [D, warpedAux] = diffeomorphicDemon(source, target, aux)
% diffeomorphicDemon applies diffeomorphic demons registration
% from 'source' to 'target', and warps 'aux' with this transformation.

% args:
% source: moving image (double, grayscale in [0,1])
% target: fixed image (double, grayscale in [0,1])
% aux: image to which apply the transformation

% returns:
% D: displacement field from source to target
% warpedAux: aux warped into target
if size(source, 3) == 3, source = rgb2gray(source); end
if size(target, 3) == 3, target = rgb2gray(target); end
source = im2double(source);
target = im2double(target);

% Diffeomorphic demons
iters = [5];
accSmooth = 2.0;

[D, ~] = imregdemons(source, target, iters, ...
    "AccumulatedFieldSmoothing", accSmooth, ...
    "PyramidLevels", numel(iters), ...
    "DisplayWaitbar", false);

warpedAux = imwarp(aux, D, "nearest");

%freeze pixel where warp is minimal
mask = hypot(D(:, :, 1), D(:, :, 2)) < 0.5;
warpedAux(mask) = source(mask);
end

function saveStabilizedVideoGif(arr, filename)

for k = 1:size(arr, 3)
    frame = arr(:, :, k);
    frame = mat2gray(frame);
    [A, map] = gray2ind(frame, 256);

    if k == 1
        imwrite(A, map, filename, "gif", "LoopCount", Inf, "DelayTime", 0.05);
    else
        imwrite(A, map, filename, "gif", "WriteMode", "append", "DelayTime", 0.05);
    end

end

end

function saveAngleAndNormOfDisplacementField(arr, angleFilename, normFilename)

for k = 1:size(arr, 4)
    D = arr(:, :, :, k);
    Dc = complex(D(:, :, 1), D(:, :, 2));

    %map abs(Dc) to the min and max of the value of the image
    %maybe comapre to min and max of the video not the image to keep quantitative value idk
    absLog = log1p(abs(Dc));
    absLog = mat2gray(absLog);
    absLog = imadjust(absLog, [prctile(absLog(:), 1) prctile(absLog(:), 99); ], [0 1]);

    % Convert to indexed
    [AAbs, mapAbs] = gray2ind(absLog, 256);
    [AAngle, mapAngle] = gray2ind(angle(Dc), 256);

    % Write to GIF
    if k == 1
        imwrite(AAbs, mapAbs, normFilename, "gif", "LoopCount", Inf, "DelayTime", 0.05);
        imwrite(AAngle, mapAngle, angleFilename, "gif", "LoopCount", Inf, "DelayTime", 0.05);
    else
        imwrite(AAbs, mapAbs, normFilename, "gif", "WriteMode", "append", "DelayTime", 0.05);
        imwrite(AAngle, mapAngle, angleFilename, "gif", "WriteMode", "append", "DelayTime", 0.05);
    end

end

end

function I = safeConvertFrame(frame)

if size(frame, 3) == 3
    I = rgb2gray(frame);
else
    I = frame;
end

I = im2double(I); % keep everything in [0,1]
end
