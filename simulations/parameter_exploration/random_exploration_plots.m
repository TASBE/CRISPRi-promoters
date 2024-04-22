% random_exploration_plots.m
% Last updated: 2023-08-21, by Helen Scott
%
% REVISION HISTORY:
%   2023-09-17 - Helen Scott
%       * Revert to run all models
%   2023-08-21 - Helen Scott
%       * Change code to only run on the multisite (4 sites) model
%   12/06/2021 - Helen Scott
%       * Copied over fuxx_perturbation_plots.m from NIH-CRISPR and made
%           necessary changes to get it to run
%

% Load the model catalog
model_catalog;

% Load the parameters
base_parameters;
base_interference_matrix;

parameters = unpackIntMatrix(parameters, i_matrix);

% Get the names of all the parameters because we want to perturb every 
% parameter individually
parameterNames = keys(parameters);

% Set color bar limits for all plots
cLimits = [0 500];

% Set number of bins for all plots
nBins = [200 50];

% Plot just the colorbar
figure('visible', 'off', 'PaperUnits','inches', 'PaperPosition', [0 0 6 1.5]); % Make but don't show the figure
cb = colorbar('southoutside');
caxis(cLimits);
ylabel(cb,'Number of Trajectories Per Cell');
saveas(gcf, './random-perturbation-results/colorbar.png', 'png');

% Loop through all the models
for i=1:n_models
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};

    % Print update on model to screen
    disp(['Currently on: ', modelName]);


    % Clean up name (for making directory/files)
    cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
        ' ', '_'), '/', '-');

    % Path to the results file
    resultsPath = ['./random-perturbation-results/', cleanModelName, '/'];

    % Set the output path for all of the plots
    outpath = [resultsPath, 'plots/'];

    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end
    
    % Load the perturbation results
    load([resultsPath, 'random-perturb-all.mat'])
        
    % Plot
    densityheatmap(x, ys, nBins, cLimits, modelName, outpath);
end