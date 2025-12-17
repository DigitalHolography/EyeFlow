function VideoNonRigidRegistering(obj, apply)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

% This function performs non-rigid registration on the video frames of M0
v = obj.M0;
[numX, numY, numFrames] = size(v);
ref_img = mean(v, 3);
low_freq = imgaussfilt(ref_img, 100);
ref_img = ref_img ./ low_freq;
ref_img = (ref_img - min(ref_img(:))) / (max(ref_img(:)) - min(ref_img(:)));

stabilized = zeros(numX, numY, numFrames);
field = zeros(numX, numY, 2, numFrames);
diff = zeros(numX, numY, numFrames);

%smoothVideo = imgaussfilt3(obj.M0, [0.1 0.1 2]);
parfor k = 1:numFrames
    %get the frame, stabilize it, save it
    tgt = safeConvertFrame(v(:, :, k));
    % [f, s] = diffeoDemon(ref_img, tgt, ...
    % "STEP_SIZE", 4, ...
    % "BLOCK_SIZE", 30, ...
    % "SEARCH_RADIUS", 5, ...
    % "SSD_THRESHOLD", 0.1, ...
    % "NUM_ITERS", 1, ...
    % "LOG_LEVEL", "WARNING" ...
    % );
    [f, s] = diffeomorphicDemon(tgt, ref_img, tgt);

    field(:, :, :, k) = f;
    stabilized(:, :, k) = s;
    diff(:, :, k) = tgt - s;
end

D.diff = diff;
D.stabilized = stabilized;
D.field = field;

A1 = field(:, :, 1, :);
A2 = field(:, :, 2, :);

% Compute the temporal derivative of mag and phase
AA(:, :, :) = complex(A1, A2);
absAA = abs(AA);
[~, ~, Ft] = gradient(absAA);
D.mag_temporal_derivative = Ft;
angleAA = angle(AA);
[~, ~, Ft] = gradient(angleAA);
D.phase_temporal_derivative = Ft;

% dA1 = diff(squeeze(A1), 1, 3);
% dA2 = diff(squeeze(A2), 1, 3);
% dA = zeros(numX, numY, 2, numFrames - 1);
% dA(:, :, 1, :) = dA1;
% dA(:, :, 2, :) = dA2;

if (exportVideos)
    saveAsGifs(D);
end

obj.displacementField = D;

if apply
    % M0_ff = obj.M0_ff;
    f_RMS = obj.f_RMS;
    f_AVG = obj.f_AVG;

    dfield = D.field;

    parfor ff = 1:numFrames
        % M0_ff(:,:,ff) = imwarp(M0_ff(:,:,ff), dfield(:,:,:,ff), "nearest");
        f_RMS(:,:,ff) = imwarp(f_RMS(:,:,ff), dfield(:,:,:,ff), "nearest");
        f_AVG(:,:,ff) = imwarp(f_AVG(:,:,ff), dfield(:,:,:,ff), "nearest");
    end
    % obj.M0_ff = M0_ff;
    obj.f_RMS = f_RMS;
    obj.f_AVG = f_AVG;

end
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
iters = 3;
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

I = im2double(I);
I = (I - min(I(:))) / (max(I(:)) - min(I(:))); % keep everything in [0,1]
end

function saveAsGifs(D)
% Normalize stabilized frames to [0, 255] uint8
normalizedStabilized = zeros(size(D.stabilized), 'uint8');

for k = 1:size(D.stabilized, 3)
    frame = D.stabilized(:, :, k);
    minVal = double(min(frame(:)));
    maxVal = double(max(frame(:)));

    if maxVal == minVal
        normalizedFrame = zeros(size(frame));
    else
        normalizedFrame = (double(frame) - minVal) / (maxVal - minVal);
    end

    normalizedStabilized(:, :, k) = uint8(normalizedFrame * 255);
end

writeGifOnDisc(normalizedStabilized, 'M0_stabilized');

% Compute magnitude and phase of displacement field
Xcomp = squeeze(D.field(:, :, 1, :));
Ycomp = squeeze(D.field(:, :, 2, :));

magnitude = sqrt(Xcomp .^ 2 + Ycomp .^ 2);
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
writeGifOnDisc(mat2gray(D.mag_temporal_derivative), 'displacement_mag_temporal_derivative');
writeGifOnDisc(mat2gray(D.phase_temporal_derivative), 'displacement_phase_temporal_derivative');
writeGifOnDisc(D.diff, 'diff_M0_stab');
end
