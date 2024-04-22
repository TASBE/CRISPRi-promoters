% random_exploration
% Last updated: 2023-08-21, by Helen Scott
%
% REVISION HISTORY:
%   2023-09-17 - Helen Scott
%       * Change code to run all models together
%   2023-08-21 - Helen Scott
%       * Change code to only run on the multisite (4 sites) model
%   2023-01-31 - Helen Scott
%       * Initial implementation
%       * Goal: One push to generate results of 10000 random perturbations 
%           of tweaked parameters for all of the models
%       * Copy the "fuzz_perturbations" file from NIH-CRISPR
%       * Will use the densityheatmap.m function from NIH-CRISPR for 
%           plotting

% Load the model catalog
model_catalog;

% Load the parameters
base_parameters;
base_interference_matrix;

parameters = unpackIntMatrix(parameters, i_matrix);

% Get the names of all the parameters because we want to perturb every 
% parameter individually
parameterNames = keys(parameters);

% Loop through all the models
for i=1:n_models
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};

    % Print update on model to screen
    disp(['Currently on: ', modelName]);

    % Clean up name (for making directory/files)
    cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
        ' ', '_'), '/', '-');

    % Set the initial
    initial = containers.Map();
    initial('V1') = 10;
    initial('V2') = 3;

    % Set the timespan
    tspan = [0 100];

    % Set the number of runs to do
    nRuns = 10000; 
    
    % Set the output path
    resultsPath = './random-perturbation-results/';

    % Make an output folder for this model
    outpath = [resultsPath, cleanModelName, '/'];

    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end

    % Do the perturbation
    results = logNormalPerturbation(modelFun, parameters, ...
        parameterNames, nRuns, tspan, initial, true, false);

    % Extract just what is needed for plotting
    percentiles = results{:, 1};
    perturbedValues = results{:, 2};
    x = results{1, 4}{1};
    ys = cell(length(results), 1);
    for k = 1:length(results)
        ys{k} = results{k, 4}{3};
    end

    % Save results as a .mat file
    save([outpath, 'random-perturb-all.mat'], 'x', 'ys')
    
end
