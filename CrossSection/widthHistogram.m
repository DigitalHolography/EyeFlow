function [medianWidth, avgWidth, stdWidth] = widthHistogram(width, width_std, area, name)

ToolBox = getGlobalToolBox;

isVal = cellfun(@(x) ~isempty(x) && ~(isnumeric(x) && isnan(x)), area);
numValid = sum(isVal, 'all');

[numCircles, numBranches] = size(area);
area_mat = nan(numCircles, numBranches);

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(area{cIdx, bIdx})
            area_mat(cIdx, bIdx) = area{cIdx, bIdx};
        end

    end

end

area_mat = reshape(area_mat, 1, []);
widths = 2 * sqrt(area_mat / pi) * 1000;

% Remove outliers (beyond ±3σ)
temp_avg = mean(widths, 'omitnan');
temp_std = std(widths, 'omitnan');
valid_idx = (widths >= (temp_avg - 3*temp_std)) & (widths <= (temp_avg + 3*temp_std));
widths = widths(valid_idx);

figure("Visible", "off")
histogram(widths, 20, FaceColor = 'k', Normalization = 'probability');
hold on

medianWidth = median(widths, "omitnan");
avgWidth = mean(widths, "omitnan");
stdWidth = std(widths, "omitnan");

% Create Gaussian distribution overlay
x = linspace(0, 200, 1000);
gaussian = normpdf(x, avgWidth, stdWidth);
% Scale Gaussian to match histogram probability
gaussian = gaussian * (max(ylim)/max(gaussian)) * 0.8; 
plot(x, gaussian, 'k-', 'LineWidth', 2);

xline(medianWidth, '--', sprintf('%.0f µm', medianWidth), 'Linewidth', 2)
set(gca, 'Linewidth', 2)
pbaspect([1.618 1 1]);
xlabel("lumen cross section width (µm)")
ylabel("probability")
xlim([0 200]) % Set x-axis limits as requested

% Add annotation with μ and σ values
annotationText = sprintf('μ = %.1f µm\nσ = %.1f µm', avgWidth, stdWidth);
annotation('textbox', [0.15 0.7 0.1 0.1], 'String', annotationText, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
           'EdgeColor', 'black', 'LineWidth', 1, 'FontSize', 10);

aa = axis;
aa(4) = aa(4) * 1.14;
axis(aa);

exportgraphics(gca, fullfile(ToolBox.path_png, 'local', sprintf("%s_%s", ToolBox.folder_name, sprintf('histogram_of_%s_section_width.png', name))))

%csv output of the widths
T = table();

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(width{cIdx, bIdx})
            T.(sprintf('Width_R%d_S%d_%s', cIdx, bIdx, name)) = width{cIdx, bIdx};
            T.(sprintf('STD_Width_R%d_S%d_%s', cIdx, bIdx, name)) = width_std{cIdx, bIdx};
        end

    end

end

writetable(T, fullfile(ToolBox.path_txt, strcat(ToolBox.folder_name, '_', 'WidthTable', '_', name, '.csv')));

% New
if contains(name, 'Artery')
    ToolBox.Outputs.add('ArterialDiameterAverage', avgWidth, 'µm');
    ToolBox.Outputs.add('ArterialDiameterMedian', medianWidth, 'µm');
    ToolBox.Outputs.add('ArterialDiameterSpread', stdWidth, 'µm');
    ToolBox.Outputs.add('ArterialValidSections', numValid, '');
else
    ToolBox.Outputs.add('VenousDiameterAverage', avgWidth, 'µm');
    ToolBox.Outputs.add('VenousDiameterMedian', medianWidth, 'µm');
    ToolBox.Outputs.add('VenousValidSections', numValid, '');
end

end
