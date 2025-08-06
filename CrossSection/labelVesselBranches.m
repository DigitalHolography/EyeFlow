function [labeledVessels, n] = labelVesselBranches(vesselMask, maskSection, xy_barycenter)

% Get parameters from global toolbox
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

[numX, numY] = size(vesselMask);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
numCircles = params.json.CrossSectionsAnalysis.NumberOfCircles;
dr = (r2 - r1) / numCircles;

% Create skeleton and remove center circle
skel = bwskel(vesselMask);
cercleMask = diskMask(numX, numY, r1, 'center', [x_c / numX y_c / numY]);
skel = skel & ~cercleMask;

% Remove branch points from skeleton
branchPoints = bwmorph(skel, 'branchpoints');
skelNoBranchesPoints = skel & ~imdilate(branchPoints, strel('disk', 2));

% Remove small branches
skelNoBranchesPoints = bwareaopen(skelNoBranchesPoints, 10); % Remove small branches
[label, n] = bwlabel(skelNoBranchesPoints); % Label disconnected branchess

% Compute distance transform (vessel thickness)
D = -bwdist(~vesselMask); % Negative transform for watershed
D(~(vesselMask & maskSection)) = -Inf; % Force background to -Inf

% Use skeleton branches as markers
markers = label > 0;
markers = imdilate(markers, strel('disk', 1)); % Ensure markers are connected

% Apply watershed
L = double(watershed(imimposemin(D, markers)));
L = L .* vesselMask; % Keep only vessel regions

labeledVessels = zeros(size(vesselMask));

for i = 1:n
    branchPixels = (L == i);
    labeledVessels(branchPixels) = i;
end

labeledVessels = labeledVessels .* maskSection;
[labeledVessels, ~] = bwlabel(labeledVessels); % Final labeling

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
