function [radiusQ, radiusQSE] = averageRadiusFlux(Qcell, QSEcell)
% AVERAGERADIUSFLUX Calculate average flux and standard error across radii
%   [radiusQ, radiusQSE] = averageRadiusFlux(Qcell, QSEcell) computes the
%   average flux and standard error across radii (circles) from cell arrays
%   containing flux measurements from multiple branches at each radius.
%
% Inputs:
%   Qcell    - Cell array of flux measurements (numCircles × numBranches)
%   QSEcell  - Cell array of flux standard errors (numCircles × numBranches)
%
% Outputs:
%   radiusQ   - Average flux per radius (numCircles × numFrames)
%   radiusQSE - Resulting standard error per radius (numCircles × numFrames)

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
radiusQ = zeros(numCircles, maxFrames);
radiusQSE = zeros(numCircles, maxFrames);

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(Qcell{cIdx, bIdx})
            currentFrames = length(Qcell{cIdx, bIdx});
            framesToUse = 1:currentFrames;

            % Handle cases where current data is shorter than maxFrames
            if currentFrames < maxFrames
                radiusQ(cIdx, framesToUse) = radiusQ(cIdx, framesToUse) + Qcell{cIdx, bIdx};
                radiusQSE(cIdx, framesToUse) = radiusQSE(cIdx, framesToUse) + (QSEcell{cIdx, bIdx} .^ 2);
            else
                radiusQ(cIdx, :) = radiusQ(cIdx, :) + Qcell{cIdx, bIdx};
                radiusQSE(cIdx, :) = radiusQSE(cIdx, :) + (QSEcell{cIdx, bIdx} .^ 2);
            end

        end

    end

end

radiusQSE = sqrt(radiusQSE);

end
