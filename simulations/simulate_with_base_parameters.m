% This file generates paper compares of all circuits with basal parameters

base_parameters;
base_interference_matrix;
model_catalog;

parameters = unpackIntMatrix(parameters, i_matrix);

initial = containers.Map();
initial('V1') = 10;
initial('V2') = 3;

final_time = 200;
time = [0 final_time];
species_names = cell(n_models, 1);
y_out = nan(n_models, final_time+1);
y_complete = cell(n_models,1);
fprintf('Simulating with base parameters');
for i=1:n_models
    try
        [time_interval, species_names{i}, y_out(i,:), y_complete{i}] = models{i,MODEL_FUN}(time,parameters,initial,1);
        fprintf('%s. \n', models{i});
    catch
        fprintf('%s! \n', models{i});
    end
end

%% Plot GFP vs time
% Save the different line colors for easy re-use
lightGreen = [0.8242 0.8398 0.6953];
darkGreen = [0.2706 0.5098 0.4314];
gold = [0.9843 0.7176 0.0510];
orange = [0.8711 0.5078 0.2656];

% Line types and colors are hard-coded to work with 12 models (one base
% exprssion, one 1 gRNA, and 2-6 gRNAs for the multiplexed and multisite
% models).
% FIXME: the "o" line type is totally overlapped, and just looks like a
% solid line
line_types = ["-", "--", ... % For the base expression and single gRNA
              "plotDash", "-.", ":", "--", "-", ... % For the multiplexed lines
              "plotDash", "-.", ":", "--", "-"]; % For the multisite lines
colors = [lightGreen; % Light green for base expression
          lightGreen; % Light green for 1 gRNA
          darkGreen; % Green for Multiplexed 2
          darkGreen; % Green for Multiplexed 3
          darkGreen; % Green for Multiplexed 4
          darkGreen; % Green for Multiplexed 5
          darkGreen; % Green for Multiplexed 6
          gold; % Gold for Multi-site 2
          gold; % Gold for Multi-site 3
          gold; % Gold for Multi-site 4
          gold; % Gold for Multi-site 5
          gold % Gold for Multi-site 6
         ];
h = figure('PaperPosition',[1 1 5, 6.5]); 
for i=1:n_models
    xPlot = time_interval+parameters('initial_delay');
    yPlot = y_out(i,:);
    % For the one line that says "plotDash" for the lineSpec use the
    % plotDash function to plot
    lineSpec = line_types(i);
    if ~any(strcmp(["-", "--", "-.", ":"], lineSpec))
        p_0 = plotDash(h, xPlot, yPlot, [1, 0.1 0.2 0.1], 20, 'Color', colors(i,:), 'LineWidth', 1);
    else
        plot(xPlot, yPlot, lineSpec, 'Color', colors(i,:), 'LineWidth', 1);
    end
    
    hold on;
end

% Style
legend_names = {models{:, 1}};
% The legend does not plot the plotDash line correctly (expected), and I
% could not figure out how to use the plotLegendLine function to plot
% outside the plot area. I wanted to edit the legned handle object, but I
% couldn't figure out how to plot the plotDash lines in the legend, only the
% native line types. So I removed the lines, and will add them in in
% PowerPoint.
[lgdh, objh] = legend(legend_names,'Location','SouthOutside');
% I tried to retrieve the object for a specific line using:
% findobj(objh, 'Type', 'Line', 'DisplayName', '2 Heterogeneous Target Sites')
% But it only returned empty arrays, so hardcoding
heteroLineH = objh(17);
identicalLineH = objh(27);
% Remove the lines
set(heteroLineH, 'LineStyle', 'none')
set(identicalLineH, 'LineStyle', 'none')
% Other style
fontsize(objh, 6, 'points');
xlabel('Hours');
set(gca, 'YScale', 'log')
ylabel('[GFP]');

% Save
outputfig(h,'base_parameters_gfp','plots');

%% Plot mean fold repression as bars
% I want the models to show up as [1 2h 2i 3h 3i 4h 4i 5h 5i 6h 6i] which
% is model index [2 3 8 4 9 5 10 6 11 7 12]
% Calculate the mean fold change for each model between time points 72-96
% Check the values to take out based on the x value
% disp(xPlot(55:79))
mean_fold_changes = zeros(n_models/2, 2);

for i = 1:n_models - 1
    if i <= 6
        idx = i;
    elseif i > 6
        idx = i + 1;
    end
    fold_change_timecourse = y_out(1, 55:79) ./ y_out(i + 1, 55:79); % TODO: CHECK WITH JAKE!
    mean_fold_changes(idx) =  mean(fold_change_timecourse, "omitnan");
end

% Plot the mean fold repression as bars
h = figure('PaperPosition',[1 1 5 6.5]);
b = bar(mean_fold_changes, 'grouped');
hold on

% Specify Colors
b(1).FaceColor = 'flat';
% Make the first bar (1 gRNA light green)
b(1).CData(1, :) = lightGreen;
% Make the rest of the first series (Heterogeneous target sites) dark green
b(1).CData(2:end, :) = repmat(darkGreen, 5, 1);
% Make the second series (Identical target sites) gold
b(2).FaceColor = gold;

% Set y axis to be log and give it a label
ylabel('Fold Repression');
set(gca, 'YScale', 'log') % Log axis makes the first bar dissapear
ylim([1 100]) % Set the lower-limit as 1 since log10(1) is 0

% Make and give labels by the number of target sites
grnaOptions = [1, 2, 3, 4, 5, 6];
xLabs = cell(size(grnaOptions));
for i = 1:length(grnaOptions)
    xLabs{i} = [num2str(grnaOptions(i)) ' Target Sites'];
end
set(gca, 'xticklabel', xLabs)
xtickangle(90)

% Key for colors
legendLabels = {'1 Target Site', 'Heterogenous Target Sites', 'Identical Target Sites'};
legendColors = {lightGreen, darkGreen, gold};

legendX = 0.75;
legendYStart = 85;
yOffset = 15;

for i = 1:length(legendLabels)
    % Add the legend label
    text(legendX, legendYStart - (i-1)*yOffset, legendLabels{i}, 'FontSize', 8, 'Color', legendColors{i})
end

% Save
outputfig(h,'base_parameters_fold_repression','plots');

%% Plot all the species for each model
% Loop through models
for i=1:n_models
    % Get the names for all the species in the model
    legend_names = species_names{i};
    % Get the results for all species in the model
    full_results = y_complete{i};
    % Plot all the species (on one graph?)
    h = figure(); 
    for spIndex = 1:length(legend_names)
        plot(time_interval+parameters('initial_delay'), full_results(spIndex, :))
        hold on;
    end
    % Plot style
    title(['All Species Concentrations for ', models{i}])
    lgd = legend(legend_names,'Location','SouthOutside');
    set(gca, 'YScale', 'log')
    % Save plot
    outputfig(h,['all_species_', regexprep(models{i}, ' ', '_')],'plots');
end