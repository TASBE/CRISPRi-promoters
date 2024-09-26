% run_all_param_scans
% Run simulations of the three models (1 gRNA, and 4 gRNA/4 site) with
% evenly spaced perturbations to each parameter independently.
% Plots the "double rainbow" plot of the gRNA on/off and the "rainbow"
% plot of fold repression.

%% Set-Up
% Load the models and all the parameters
model_catalog;
base_parameters;
base_interference_matrix;
parameters = unpackIntMatrix(parameters, i_matrix);

% Set the output directory
outpath = './results/';

% Make a container map of the parameters I want to scan with the range to
% use (i.e. 1 order of magnitude lower, 2 orders of magnitude higher for
% dCas9 degradation).
parametersToScan = containers.Map();
parametersToScan('Cas_gRNA_binding') = [2 2];
parametersToScan('alpha_r_gRNA1') = [2 2];
parametersToScan('alpha_p_GFP') = [2 2];
parametersToScan('alpha_p_dCas9') = [2 2];
parametersToScan('delta_g') = [2 2];
parametersToScan('lambda') = [2 1];
parametersToScan('K_R') = [2 2];

% Settings for the model simulations
% Set the timespan
tspan = [0 200];

% Make a dictionary of the description of the V1 value (a string) as the
% key, and the V1 value to use as the value
v1Values = [0, 10];
v1names = ["Off", "On"];
v1Levels = containers.Map(v1names, v1Values);

% For each model
for i = 1:n_models
    % Print the model name to keep track of progress
    disp(strcat("Currently on: ", models{i}))

    % Add model specific variables to the parameter scan
    for j = 1:length(models{i, MODEL_PARAMS})
        parametersToScan(models{i, MODEL_PARAMS}(j)) = [2 2];
    end

    % Call function to preform scan
    scanParameters(models(i, :), v1Levels, tspan, parameters, parametersToScan, outpath)
    
    % Call function to calcuate fold repression
    calculateFoldRepression(models(i, :), parametersToScan, outpath)
    
    % Call function to make double rainbow plots
    plotMEFL(models(i, :), parameters, parametersToScan, outpath)

    % Call function to make rainbow plot
    plotFoldRepression(models(i, :), parameters, parametersToScan, outpath)
end
