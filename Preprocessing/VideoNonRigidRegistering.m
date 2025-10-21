function VideoNonRigidRegistering(obj)
% This function performs non-rigid registration on the video frames of M0
v = obj.M0;
[numX, numY, numFrames] = size(v);
ref_img = mean(v, 3);
low_freq = imgaussfilt(ref_img, 35);
ref_img = ref_img ./ low_freq;

stabilized = zeros(numX, numY, numFrames);

%smoothVideo = imgaussfilt3(obj.M0, [0.1 0.1 2]);
parfor k = 1:numFrames
    %get the frame, stabilize it, save it
    tgt = safeConvertFrame(v(:, :, k));
    tgt = tgt ./ low_freq;
    [~, s] = diffeomorphicDemon(tgt, ref_img, tgt);

    stabilized(:, :, k) = s;
end

ref_img2 = log(mean(stabilized, 3));
field = zeros(numX, numY, 2, numFrames);

parfor k = 1:numFrames
    %get the frame, stabilize it, save it
    tgt = safeConvertFrame(v(:, :, k));
    tgt = tgt ./ low_freq;
    [f, ~] = diffeomorphicDemon(tgt, ref_img2, tgt);

    field(:, :, :, k) = f;
end

D.stabilized = stabilized;
D.field = field;

A1 = field(:, :, 1, :);
A2 = field(:, :, 2, :);

AA(:,:,:) =complex(A1,A2);
absAA = abs(AA);
[~, ~, Ft] = gradient(absAA);

% dA1 = diff(squeeze(A1), 1, 3);
% dA2 = diff(squeeze(A2), 1, 3);
% dA = zeros(numX, numY, 2, numFrames - 1);
% dA(:, :, 1, :) = dA1;
% dA(:, :, 2, :) = dA2;  

D.gradient = Ft;

saveAsGifs(D);
obj.displacementField = D;
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
iters = 5;
accSmooth = 1.0;

[D, ~] = imregdemons(source, target, iters, ...
    "AccumulatedFieldSmoothing", accSmooth, ...
    "PyramidLevels", numel(iters), ...
    "DisplayWaitbar", false);

warpedAux = imwarp(aux, D, "nearest");

%freeze pixel where warp is minimal
mask = hypot(D(:, :, 1), D(:, :, 2)) < 0.5;
warpedAux(mask) = source(mask);
end

function I = safeConvertFrame(frame)

if size(frame, 3) == 3
    I = rgb2gray(frame);
else
    I = frame;
end

I = im2double(I); % keep everything in [0,1]
end

function saveAsGifs(D)
% Normalize stabilized frames to [0, 255] uint8
normalizedStabilized = zeros(size(D.stabilized), 'uint8');
for k = 1:size(D.stabilized, 3)
    frame = D.stabilized(:,:,k);
    minVal = double(min(frame(:)));
    maxVal = double(max(frame(:)));
    if maxVal == minVal
        normalizedFrame = zeros(size(frame));
    else
        normalizedFrame = (double(frame) - minVal) / (maxVal - minVal);
    end
    normalizedStabilized(:,:,k) = uint8(normalizedFrame * 255);
end
writeGifOnDisc(normalizedStabilized, 'stabilized');

% Compute magnitude and phase of displacement field
Xcomp = squeeze(D.field(:, :, 1, :));
Ycomp = squeeze(D.field(:, :, 2, :));

magnitude = sqrt(Xcomp.^2 + Ycomp.^2);
phase = atan2(Ycomp, Xcomp);

% Normalize magnitude to [0,255] uint8
magMin = min(magnitude(:));
magMax = max(magnitude(:));
if magMax == magMin
    normMagnitude = zeros(size(magnitude));
else
    normMagnitude = (magnitude - magMin) / (magMax - magMin);
end
normMagnitude = uint8(normMagnitude * 255);
writeGifOnDisc(normMagnitude, 'displacement_magnitude');

% Normalize phase from [-pi, pi] to [0, 255] uint8
normPhase = (phase + pi) / (2 * pi); % normalize to [0,1]
normPhase = uint8(normPhase * 255);
writeGifOnDisc(normPhase, 'displacement_phase');
writeGifOnDisc(mat2gray(D.gradient), 'displacement_temporal_derivative');
end