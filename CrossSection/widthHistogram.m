function widthHistogram(width, width_std, area, name)

ToolBox = getGlobalToolBox;

figure("Visible", "off")
histogram(2 * sqrt(area(area ~= 0) / pi) * 1000, 50, FaceColor = 'k');
hold on

medianWidth = median(2 * sqrt(area(area ~= 0) / pi) * 1000, "omitnan");
avgWidth = mean(2 * sqrt(area(area ~= 0) / pi) * 1000, "omitnan");
stdWidth = std(2 * sqrt(area(area ~= 0) / pi) * 1000, "omitnan");
xline(medianWidth, '--', sprintf('%.0f µm', medianWidth), 'Linewidth', 2)
set(gca, 'Linewidth', 2)

aa = axis;
aa(4) = aa(4) * 1.14;
axis(aa);
title(sprintf('Histogram of %s sections width (µm)', name));

exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_%s", ToolBox.main_foldername, sprintf('histogram_of_%s_section_width.png', name))))

%csv output of the widths
T = table();
numR = length(width); % number of radii

for rIdx = 1:numR
    numSection = length(width{rIdx});

    for sectionIdx = 1:numSection
        T.(sprintf('Width_R%d_S%d_%s', rIdx, sectionIdx, name)) = squeeze(squeeze(width{rIdx}(sectionIdx)));
        T.(sprintf('STD_Width_R%d_S%d_%s', rIdx, sectionIdx, name)) = squeeze(squeeze(width_std{rIdx}(sectionIdx)));
    end

end

writetable(T, fullfile(ToolBox.path_txt, strcat(ToolBox.main_foldername, '_', 'WidthTable', '_', name, '.csv')));


% New
if contains(name, 'Artery')
    ToolBox.Outputs.add('ArterialDiameterAverage') = avgWidth;
    ToolBox.Outputs.add('ArterialDiameterSpread') = stdWidth;
else
    ToolBox.Outputs.add('VenousDiameterAverage') = avgWidth;
    ToolBox.Outputs.add('VenousDiameterSpread') = stdWidth;
end

end
