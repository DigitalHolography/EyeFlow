function [radius_Q, radius_Q_SE] = averageRadiusFlux(Q_cell, Q_SE_cell)
% AVERAGERADIUSFLUX Calculate average flux and standard error across radii
%   [radius_Q, radius_Q_SE] = averageRadiusFlux(Q_cell, Q_SE_cell) computes the
%   average flux and standard error across radii (circles) from cell arrays
%   containing flux measurements from multiple branches at each radius.
%
% Inputs:
%   Q_cell    - Cell array of flux measurements (numCircles × numBranches)
%   Q_SE_cell  - Cell array of flux standard errors (numCircles × numBranches)
%
% Outputs:
%   radius_Q   - Average flux per radius (numCircles × numFrames)
%   radius_Q_SE - Resulting standard error per radius (numCircles × numFrames)

% Validate input dimensions
[numCircles, numBranches] = size(Q_cell);

if ~isequal(size(Q_SE_cell), [numCircles, numBranches])
    error('Q_cell and Q_SE_cell must have the same dimensions');
end

% Determine number of frames (maximum length in all cells)
allLengths = cellfun(@length, Q_cell(:));
numFrames = max(allLengths(:));

% Initialize outputs
radius_Q = zeros(numCircles, numFrames);
radius_Q_SE = zeros(numCircles, numFrames);

% Compute sums for each radius
for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(Q_cell{cIdx, bIdx})

            radius_Q(cIdx, :) = radius_Q(cIdx, :) + Q_cell{cIdx, bIdx};
            radius_Q_SE(cIdx, :) = radius_Q_SE(cIdx, :) + (Q_SE_cell{cIdx, bIdx} .^ 2);

        end

    end

end

radius_Q_SE = sqrt(radius_Q_SE);

end
