%% Load the results from a level where things are strongly different (e.g.10^4 alpha_p_dCas9)

% Load the results from the multiplexed-gRNA model
load("parameter_exploration/parameter-scan-results/multiplexed_grna_repression/no_interaction/scan-alpha_p_dCas9-gRNAs-On-On.mat")
multiplexed_results = results{end, 3};

% Load the results from the multi-site single-gRNA model
load("parameter_exploration/parameter-scan-results/multisite_repression/one_gRNA/scan-alpha_p_dCas9-gRNAs-On.mat")
multisite_results = results{end, 3};

%% Plot all of the concentrations (not just gRNA)
figure;
hold on;

% Have the hardcode the names of each variable (they are different for the
% different models) and are not save in the results file
multiplexed_ys = {'GFP', 'V1', 'V2', 'dCas9', 'dCas9_gRNA1', 'dCas9_gRNA2', 'gRNA1', 'gRNA2'};

% Make color map of 8 colors
cmap = hsv(8);

% Plot all the concentrations in the multiplexed results
for i = 1:height(multiplexed_results{3})
    label = strcat(multiplexed_ys{i}, ' (Multiplexed)');
    concentration_list = multiplexed_results{3}(i, :);
    plot(multiplexed_results{1}, concentration_list, ...
        'Color', cmap(i, :), 'LineStyle', '-', ...
        'DisplayName', label)
end

% Add a line for the total concentration of dCas9 (unbound, and both
% complexes)
plot(multiplexed_results{1}, sum(multiplexed_results{1,3}(4:6, :)), ...
        'Color', 'black', 'LineStyle', '-', ...
        'DisplayName', 'Total dCas9 (Multiplexed)')

multisite_ys = {'GFP', 'V1', 'V2', 'dCas9', 'dCas9_gRNA1', 'gRNA1'};

% Plot all the concentrations in the multisite results
for i = 1:height(multisite_results{3})
    label = strcat(multisite_ys{i}, ' (Multi-site)');
    color_value = cmap(find(strcmp(multiplexed_ys, multisite_ys{i})), :);
    concentration_list = multisite_results{3}(i, :);
    plot(multisite_results{1}, concentration_list, ...
        'Color', cmap(i, :), 'LineStyle', '--', ...
        'DisplayName', label)
end

% Add a line for the total concentration of dCas9 (unbound, and both
% complexes)
plot(multiplexed_results{1}, sum(multisite_results{1,3}(4:5, :)), ...
        'Color', 'black', 'LineStyle', '--', ...
        'DisplayName', 'Total dCas9 (Multisite)')

% Style
legend();
xlabel('Time (hours)');
ylabel('Concentration');
set(gca, 'YScale', 'log')


% Save the figure with the default bounds
saveas(gcf, 'compare_all_concentrations.png')