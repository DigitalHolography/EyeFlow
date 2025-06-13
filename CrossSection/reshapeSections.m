function [A_mat, Q_mat, Q_std_mat, Q_branch, Q_se_branch] = reshapeSections(numFrames, numSections, A_cell, Q_cell, Q_se_cell)

numSectionsMax = max(numSections);
[numCircles, numBranches] = size(numSections, 2);

A_mat = zeros(numCircles, numSectionsMax);
Q_mat = zeros(numCircles, numSectionsMax, numFrames);
Q_std_mat = zeros(numCircles, numSectionsMax, numFrames);

for cIdx = 1:numCircles
    numSection = numSections(cIdx);

    if (numSection < numSectionsMax) && (numSection ~= 0)

        for sectionIdx = numSection:numSectionsMax
            A_mat(cIdx, sectionIdx) = nan;
        end

    end

    if numSection ~= 0

        for sectionIdx = 1:numSection
            A_mat(cIdx, sectionIdx) = A_cell{cIdx}(sectionIdx);
        end

    end

end

end
