function [radius_v, radius_v_SE] = averageRadiusVelocity(v_cell, v_SE_cell)
% AVERAGERADIUSVELOCITY Calculate average velocity and standard error across radii
%   [radius_v, radius_v_SE] = averageRadiusVelocity(v_cell, v_SE_cell) computes the
%   average velocity and standard error across radii (circles) from cell arrays
%   containing velocity measurements from multiple branches at each radius.
%
% Inputs:
%   v_cell    - Cell array of velocity measurements (numCircles × numBranches)
%   v_SE_cell  - Cell array of velocity standard errors (numCircles × numBranches)
%
% Outputs:
%   radius_v   - Average velocity per radius (numCircles × numFrames)
%   radius_v_SE - Resulting standard error per radius (numCircles × numFrames)

% Validate input dimensions
[numCircles, numBranches] = size(v_cell);

if ~isequal(size(v_SE_cell), [numCircles, numBranches])
    error('v_cell and v_SE_cell must have the same dimensions');
end

% Determine number of frames (maximum length in all cells)
allLengths = cellfun(@length, v_cell(:));
numFrames = max(allLengths(:));

% Initialize outputs
radius_v = zeros(numCircles, numFrames);
radius_v_SE = zeros(numCircles, numFrames);

% Compute sums for each radius
for cIdx = 1:numCircles

    numBranches_per_Circle = 0; % Track the number of branches per circle

    for bIdx = 1:numBranches

        if ~isempty(v_cell{cIdx, bIdx})

            radius_v(cIdx, :) = radius_v(cIdx, :) + v_cell{cIdx, bIdx};
            radius_v_SE(cIdx, :) = radius_v_SE(cIdx, :) + (v_SE_cell{cIdx, bIdx} .^ 2);
            numBranches_per_Circle = numBranches_per_Circle + 1;

        end

    end

    % Average over branches
    if numBranches_per_Circle > 0
        radius_v(cIdx, :) = radius_v(cIdx, :) ./ numBranches_per_Circle;
        radius_v_SE(cIdx, :) = sqrt(radius_v_SE(cIdx, :)) ./ numBranches_per_Circle;
    end

end

end
