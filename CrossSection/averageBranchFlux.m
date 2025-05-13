function [branchQ, branchQSE] = averageBranchFlux(Qcell, QSEcell)
% AVERAGEBRANCHFLUX Calculate average flux and standard error across branches
%   [branchQ, branchQSE] = averageBranchFlux(Qcell, QSEcell) computes the
%   average flux and standard error across branches from cell arrays containing
%   flux measurements from multiple circles.
%
% Inputs:
%   Qcell    - Cell array of flux measurements (numCircles × numBranches)
%   QSEcell  - Cell array of flux standard errors (numCircles × numBranches)
%
% Outputs:
%   branchQ   - Average flux per branch (numBranches × numFrames)
%   branchQSE - Resulting standard error per branch (numBranches × numFrames)

[numCircles, numBranches] = size(Qcell);

if ~isequal(size(QSEcell), [numCircles, numBranches])
    error('Qcell and QSEcell must have the same dimensions');
end

% Determine maximum number of frames across all data
maxFrames = 0;

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(Qcell{cIdx, bIdx})
            maxFrames = max(maxFrames, length(Qcell{cIdx, bIdx}));
        end

    end

end

if maxFrames == 0
    error('Input cells contain no valid data');
end

% Initialize outputs
branchQ = zeros(numBranches, maxFrames);
branchQSE = zeros(numBranches, maxFrames);
validCounts = zeros(numBranches, maxFrames); % Track valid counts per frame

for bIdx = 1:numBranches

    for cIdx = 1:numCircles

        if ~isempty(Qcell{cIdx, bIdx})
            currentFrames = length(Qcell{cIdx, bIdx});

            % Pad with NaN if needed (for shorter sequences)
            if currentFrames < maxFrames
                framesToUse = 1:currentFrames;
                branchQ(bIdx, framesToUse) = branchQ(bIdx, framesToUse) + Qcell{cIdx, bIdx};
                branchQSE(bIdx, framesToUse) = branchQSE(bIdx, framesToUse) + (QSEcell{cIdx, bIdx} .^ 2);
                validCounts(bIdx, framesToUse) = validCounts(bIdx, framesToUse) + 1;
            else
                branchQ(bIdx, :) = branchQ(bIdx, :) + Qcell{cIdx, bIdx};
                branchQSE(bIdx, :) = branchQSE(bIdx, :) + (QSEcell{cIdx, bIdx} .^ 2);
                validCounts(bIdx, :) = validCounts(bIdx, :) + 1;
            end

        end

    end

end

% Calculate averages only where we have data
for bIdx = 1:numBranches
    hasData = validCounts(bIdx, :) > 0;

    if any(hasData)
        branchQ(bIdx, hasData) = branchQ(bIdx, hasData) ./ validCounts(bIdx, hasData);
        branchQSE(bIdx, hasData) = sqrt(branchQSE(bIdx, hasData)) ./ validCounts(bIdx, hasData);
    end

    % Leave as zeros where no data exists
end

end
