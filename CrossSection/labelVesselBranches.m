function [labeledVessels, n] = labelVesselBranches(vesselMask, maskSection, xy_barycenter, fun_params)
% labelVesselBranches - Label vessel branches in a binary vessel mask using watershed segmentation.
%
% Syntax:
%   [labeledVessels, n] = labelVesselBranches(vesselMask, maskSection, xy_barycenter)
%
% Inputs:
%   vesselMask     - A binary mask (logical array) of the segmented vessels.
%   maskSection    - A binary mask (logical array) defining the cross-sectional area of interest.
%   xy_barycenter  - A 2-element vector [x, y] representing the center of the cross-section.
%
% Outputs:
%   labeledVessels - A labeled mask (integer array) where each vessel branch has a unique label.
%   n              - The number of labeled vessel branches. 0 if none are found.

arguments
    vesselMask logical % Binary mask of the segmented vessels
    maskSection logical % Binary mask defining the cross-sectional area of interest
    xy_barycenter (1, 2) double % [x, y] coordinates of the center of the cross-section
    fun_params.refine logical = true % perform small refining or not and apply a sectioning
    fun_params.min_area_percent double = 0.1 % area percentage of the image minimal for each label
    fun_params.strelSize double = 3 % Size for structuring element in morphological operations
    fun_params.smallBranchCriteria double = 10 % Minimum area (in pixels) for a branch to be considered valid
end

% Get parameters from global ToolBox
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% Get image size and define parameters
[numX, numY] = size(vesselMask);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
numCircles = params.json.generateCrossSectionSignals.NumberOfCircles;
dr = (r2 - r1) / numCircles;
SE = strel('square', fun_params.strelSize);

% Create skeleton and remove center circle
skel = bwskel(vesselMask);

% Remove branch points from skeleton
branchPoints = bwmorph(skel, 'branchpoints');
skelNoBranchesPoints = skel & ~imdilate(branchPoints, SE);

% Remove small branches
skelNoBranchesPoints = bwareaopen(skelNoBranchesPoints, fun_params.smallBranchCriteria); % Remove small branches

% Label remaining branches
[label, n] = bwlabel(skelNoBranchesPoints);

% Compute distance transform (vessel thickness)
D = -bwdist(~vesselMask); % Negative transform for watershed
D(~(vesselMask & maskSection)) = -Inf; % Force background to -Inf

% Use skeleton branches as markers
markers = label > 0;
% markers = imdilate(markers, SE); % Ensure markers are connected

% Apply watershed
L = double(watershed(imimposemin(D, markers)));
L = L .* vesselMask; % Keep only vessel regions

labeledVessels = L;

% minAreaThreshold = floor(numX * numY * fun_params.min_area_percent / 100);
% labeledVessels = bwareaopen(labeledVessels, minAreaThreshold);

if ~fun_params.refine
    return;
end

cercleMask = diskMask(numX, numY, r1, r2, 'center', [x_c / numX y_c / numY]);
labeledVessels = labeledVessels .* cercleMask;
[labeledVessels, n] = bwlabel(labeledVessels); % Final labeling



% Remove small spots
labeledVesselsClean = false(numX, numY);

for cIdx = 1:numCircles
    r_in = r1 + (cIdx - 1) * dr;
    r_out = r_in + dr;
    maskSectionCircle = diskMask(numX, numY, r_in, r_out, center = [x_c / numX, y_c / numY]);

    minAreaThreshold = floor(numX / 10); % Adjust this value as needed
    labeledVesselsClean = labeledVesselsClean | bwareaopen(labeledVessels & maskSectionCircle, minAreaThreshold);
end

[labeledVessels, n] = bwlabel(labeledVesselsClean); % Relabel after removal

end
