function [D_mid, D_avg, D_std] = widthHistogram(D, dD, A, name)

ToolBox = getGlobalToolBox;

isVal = cellfun(@(x) ~isempty(x) && ~(isnumeric(x) && isnan(x)), A);
numValid = sum(isVal, 'all');

[numCircles, numBranches] = size(A);
area_mat = nan(numCircles, numBranches);

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(A{cIdx, bIdx})
            area_mat(cIdx, bIdx) = A{cIdx, bIdx};
        end

    end

end

area_mat = reshape(area_mat, 1, []);
diameters = 2 * sqrt(area_mat / pi) * 1000;

% Remove outliers (beyond ±3σ)
temp_avg = mean(diameters, 'omitnan');
temp_std = std(diameters, 'omitnan');
valid_idx = (diameters >= (temp_avg - 3*temp_std)) & (diameters <= (temp_avg + 3*temp_std));
diameters = diameters(valid_idx);

figure("Visible", "off")
histogram(diameters, 20, FaceColor = 'k', Normalization = 'probability');
hold on

D_mid = median(diameters, "omitnan");
D_avg = mean(diameters, "omitnan");
D_std = std(diameters, "omitnan");

% Create Gaussian distribution overlay
x = linspace(0, 200, 1000);
gaussian = normpdf(x, D_avg, D_std);
% Scale Gaussian to match histogram probability
gaussian = gaussian * (max(ylim)/max(gaussian)) * 0.8; 
plot(x, gaussian, 'k-', 'LineWidth', 2);

xline(D_mid, '--', sprintf('%.0f µm', D_mid), 'Linewidth', 2)
set(gca, 'Linewidth', 2)
pbaspect([1.618 1 1]);
xlabel("lumen cross section diameter (µm)")
ylabel("probability")
xlim([0 200]) % Set x-axis limits as requested

% Add annotation with μ and σ values
annotationText = sprintf('Average = %.1f µm\nSpread = %.1f µm\nMedian = %.1f µm', D_avg, D_std, D_mid);
annotation('textbox', [0.15 0.7 0.1 0.1], 'String', annotationText, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
           'EdgeColor', 'none', 'LineWidth', 1, 'FontSize', 10);

aa = axis;
aa(4) = aa(4) * 1.14;
axis(aa);

exportgraphics(gca, fullfile(ToolBox.path_png, 'local', sprintf("%s_%s", ToolBox.folder_name, sprintf('histogram_of_%s_section_diameter.png', name))))

%csv output of the widths
T = table();

for cIdx = 1:numCircles

    for bIdx = 1:numBranches

        if ~isempty(D{cIdx, bIdx})
            T.(sprintf('Width_R%d_S%d_%s', cIdx, bIdx, name)) = D{cIdx, bIdx};
            T.(sprintf('STD_Width_R%d_S%d_%s', cIdx, bIdx, name)) = dD{cIdx, bIdx};
        end

    end

end

writetable(T, fullfile(ToolBox.path_txt, strcat(ToolBox.folder_name, '_', 'WidthTable', '_', name, '.csv')));

% New
if contains(name, 'Artery')
    ToolBox.Outputs.add('ArterialDiameterAverage', D_avg, 'µm');
    ToolBox.Outputs.add('ArterialDiameterMedian', D_mid, 'µm');
    ToolBox.Outputs.add('ArterialDiameterSpread', D_std, 'µm');
    ToolBox.Outputs.add('ArterialValidSections', numValid, '');
else
    ToolBox.Outputs.add('VenousDiameterAverage', D_avg, 'µm');
    ToolBox.Outputs.add('VenousDiameterMedian', D_mid, 'µm');
    ToolBox.Outputs.add('VenousValidSections', numValid, '');
end

end
