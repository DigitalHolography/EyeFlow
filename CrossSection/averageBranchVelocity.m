function [branch_v, branch_v_SE] = averageBranchVelocity(v_cell, v_SE_cell)
% AVERAGEBRANCHVELOCITY Calculate average velocity and standard error across branches
%   [branch_v, branch_v_SE] = averageBranchVelocity(v_cell, v_SE_cell) computes the
%   average velocity and standard error across branches from cell arrays containing
%   velocity measurements from multiple circles.
%
% Inputs:
%   v_cell    - Cell array of velocity measurements (numCircles × numBranches)
%   v_SE_cell  - Cell array of velocity standard errors (numCircles × numBranches)
%
% Outputs:
%   branch_v   - Average velocity per branch (numBranches × numFrames)
%   branch_v_SE - Resulting standard error per branch (numBranches × numFrames)

% Validate input dimensions
[numCircles, numBranches] = size(v_cell);

if ~isequal(size(v_SE_cell), [numCircles, numBranches])
    error('v_cell and v_SE_cell must have the same dimensions');
end

% Determine number of frames (maximum length in all cells)
allLengths = cellfun(@length, v_cell(:));
numFrames = max(allLengths(:));

% Initialize outputs with NaN (better than zeros for missing data)
branch_v = zeros(numBranches, numFrames);
branch_v_SE = zeros(numBranches, numFrames);

% Compute sums
for bIdx = 1:numBranches

    numCircles_per_Branch = 0; % Track the number of circles per branch

    for cIdx = 1:numCircles

        if ~isempty(v_cell{cIdx, bIdx})

            branch_v(bIdx, :) = branch_v(bIdx, :) + v_cell{cIdx, bIdx};
            branch_v_SE(bIdx, :) = branch_v_SE(bIdx, :) + (v_SE_cell{cIdx, bIdx} .^ 2);
            numCircles_per_Branch = numCircles_per_Branch + 1;

        end

    end

    % Average over circles
    if numCircles_per_Branch > 0
        branch_v(bIdx, :) = branch_v(bIdx, :) ./ numCircles_per_Branch;
        branch_v_SE(bIdx, :) = sqrt(branch_v_SE(bIdx, :)) ./ numCircles_per_Branch;
    end

end

end
