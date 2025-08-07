function [branch_Q, branch_Q_SE] = averageBranchFlux(Q_cell, Q_SE_cell)
% AVERAGEBRANCHFLUX Calculate average flux and standard error across branches
%   [branch_Q, branch_Q_SE] = averageBranchFlux(Q_cell, Q_SE_cell) computes the
%   average flux and standard error across branches from cell arrays containing
%   flux measurements from multiple circles.
%
% Inputs:
%   Q_cell    - Cell array of flux measurements (numCircles × numBranches)
%   Q_SE_cell  - Cell array of flux standard errors (numCircles × numBranches)
%
% Outputs:
%   branch_Q   - Average flux per branch (numBranches × numFrames)
%   branch_Q_SE - Resulting standard error per branch (numBranches × numFrames)

% Validate input dimensions
[numCircles, numBranches] = size(Q_cell);

if ~isequal(size(Q_SE_cell), [numCircles, numBranches])
    error('Q_cell and Q_SE_cell must have the same dimensions');
end

% Determine number of frames (maximum length in all cells)
allLengths = cellfun(@length, Q_cell(:));
numFrames = max(allLengths(:));

% Initialize outputs with NaN (better than zeros for missing data)
branch_Q = zeros(numBranches, numFrames);
branch_Q_SE = zeros(numBranches, numFrames);

% Compute sums
for bIdx = 1:numBranches

    numCircles_per_Branch = 0; % Track the number of circles per branch

    for cIdx = 1:numCircles

        if ~isempty(Q_cell{cIdx, bIdx})

            branch_Q(bIdx, :) = branch_Q(bIdx, :) + Q_cell{cIdx, bIdx};
            branch_Q_SE(bIdx, :) = branch_Q_SE(bIdx, :) + (Q_SE_cell{cIdx, bIdx} .^ 2);
            numCircles_per_Branch = numCircles_per_Branch + 1;

        end

    end

    % Average over circles
    if numCircles_per_Branch > 0
        branch_Q(bIdx, :) = branch_Q(bIdx, :) ./ numCircles_per_Branch;
        branch_Q_SE(bIdx, :) = sqrt(branch_Q_SE(bIdx, :)) ./ numCircles_per_Branch;
    end

end

end
