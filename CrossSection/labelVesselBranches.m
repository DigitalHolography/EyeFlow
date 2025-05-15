function [labeledVessels, n] = labelVesselBranches(vesselMask, maskSection, xy_barycenter)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

[numX, numY] = size(vesselMask);
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
r1 = params.json.SizeOfField.SmallRadiusRatio;

skel = bwskel(vesselMask);
cercleMask = diskMask(numX, numY, r1, 'center', [x_c / numX y_c / numY]);
skel = skel & ~cercleMask;
branchPoints = bwmorph(skel, 'branchpoints');
skelNoBranches = skel & ~imdilate(branchPoints, strel('disk', 2));
[label, n] = bwlabel(skelNoBranches); % Label disconnected branchess

% Compute distance transform (vessel thickness)
D = -bwdist(~vesselMask); % Negative transform for watershed
D(~vesselMask) = -Inf; % Force background to -Inf

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
[labeledVessels, n] = bwlabel(labeledVessels); % Label disconnected branchess

end
