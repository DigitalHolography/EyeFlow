function womersleyResultsAnalysis(womersley_results)
    arguments
        womersley_results struct % should put an array of struct later
    end

    metrics = get_array_field(womersley_results, "metrics");

    alpha_n_arr = get_array_field(womersley_results, "alpha_n");

    R0_arr = get_array_field(womersley_results, "R0");


    create_violin([alpha_n_arr]);
end


% +=====================================================================+ %
% |                          HELPER FUNCITONS                           | %
% +=====================================================================+ %

function create_violin(data, x_labels)
    arguments
        data
        x_labels string = []
    end

    figure;
    v = violinplot(data);
    hold on;

    if iscell(data)
        numGroups = length(data);
    else
        numGroups = size(data, 2);
    end

    statsText = {'\bf  Statistical Summary  '};

    for i = 1:numGroups

        if iscell(data)
            currentData = data{i};
        else
            currentData = data(:, i);
        end

        mu = mean(currentData);
        med = median(currentData);
        sigma = std(currentData);
        
        eb = errorbar(i, mu, sigma, 'r', 'LineWidth', 1.5, 'CapSize', 10);
        hm = plot(i, mu, 'rx', 'MarkerSize', 12, 'LineWidth', 3);
        hmed = plot(i, med, 'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

        % Save handles for the first iteration only (for the legend)
        if i == 1
            h_mean = hm;
            h_median = hmed;
            h_std = eb;
        end
    
        rowStr = sprintf('\\bf\\rm \\mu=%.2f, Med=%.2f, \\sigma=%.2f', ...
                         mu, med, sigma);
        statsText{end+1} = rowStr;

    end
    
    dim = [0.15 0.6 0.3 0.3]; 
    annotation('textbox', dim, 'String', statsText, ...
               'FitBoxToText', 'on', ...
               'BackgroundColor', 'white', ...
               'EdgeColor', 'black', ...
               'FaceAlpha', 0.8, ... % Slight transparency
               'FontSize', 9);

    legend([h_mean, h_median, h_std], ...
           {'Mean', 'Median', 'Standard Deviation'}, ...
           'Location', 'best');

    xticklabels(x_labels);
    ylabel('Values');
    hold off;
end

function res_array = get_array_field(womersley_results, field)
    arguments
        womersley_results struct
        field string
    end
    res_array = [womersley_results.(field)]';
    % field_parts = strsplit(field, '.');
    % 
    % current_data = womersley_results;
    % 
    % for i = 1:length(field_parts)
    %     this_field = field_parts{i};
    % 
    %     if isempty(current_data)
    %         res_array = [];
    %         return;
    %     end
    %
    %     extracted_values = {current_data.(this_field)};
    % 
    %     if i < length(field_parts)
    %         try
    %             current_data = [extracted_values{:}];
    %         catch
    %             error('Structure inconsistency found at field "%s".', this_field);
    %         end
    %     else
    %         res_array = extracted_values;
    %     end
    % end
end