%% Load the results
% Load the results from the multiplexed-gRNA model
load("vector_concentration_test_results/scan-vector-concentrations-multiplexed_gRNA_epression.mat")
multiplexed_results = results;

% Load the results from the multi-site single-gRNA model
load("vector_concentration_test_results/scan-vector-concentrations-multisite_gRNA_repression.mat")
multisite_results = results;


%% Plot
% Have the hardcode the names of each variable (they are different for the
% different models) and are not save in the results file
multiplexed_ys = {'GFP', 'V1', 'V2', 'dCas9', 'dCas9_gRNA1', 'dCas9_gRNA2', 'gRNA1', 'gRNA2'};
multisite_ys = {'GFP', 'V1', 'V2', 'dCas9', 'dCas9_gRNA1', 'gRNA1'};

% Make color map of 8 colors
cmap = hsv(8);

% For each vector concetration make a plot with the multisite and
% multiplexed concentrations
for conc_idx = 1:length(multiplexed_results)
    figure;
    hold on;

    for i = 1:height(multiplexed_results{i, 2}{3})
        label = strcat(multiplexed_ys{i}, ' (Multiplexed)');
        concentration_list = multiplexed_results{i, 2}{3}(i, :);
        plot(multiplexed_results{1}, concentration_list, ...
            'Color', cmap(i, :), 'LineStyle', '-', ...
            'DisplayName', label)
    end

    % Add a line for the total concentration of dCas9 (unbound, and both
    % complexes)
    plot(multiplexed_results{1}, sum(multiplexed_results{i, 2}{3}(4:6, :)), ...
            'Color', 'black', 'LineStyle', '-', ...
            'DisplayName', 'Total dCas9 (Multiplexed)')

    % Plot all the concentrations in the multisite results
    for i = 1:height(multisite_results{i, 2}{3})
        label = strcat(multisite_ys{i}, ' (Multi-site)');
        color_value = cmap(find(strcmp(multiplexed_ys, multisite_ys{i})), :);
        concentration_list = multisite_results{i, 2}{3}(i, :);
        plot(multisite_results{1}, concentration_list, ...
            'Color', cmap(i, :), 'LineStyle', '--', ...
            'DisplayName', label)
    end
    
    % Add a line for the total concentration of dCas9 (unbound, and both
    % complexes)
    plot(multiplexed_results{1}, sum(multisite_results{i, 2}{3}(4:5, :)), ...
            'Color', 'black', 'LineStyle', '--', ...
            'DisplayName', 'Total dCas9 (Multisite)')

    % Style
    legend();
    xlabel('Time (hours)');
    ylabel('Concentration');
    set(gca, 'YScale', 'log')
    
    
    % Save the figure with the default bounds
    saveas(gcf, ['compare_all_concentrations-vector-' conc_idx '.png'])
end