% gRNA_on_off
% Last updated: 2023-02-09, by Helen Scott
%
% REVISION HISTORY:
%   2023-02-09 - Helen Scott
%       * Renamed the file to gRNA_on_off
%       * Change from perturbing the alpha_rs to perturbing V1
%   2023-02-07 - Helen Scott
%       * Revised to run high and low alphas for just the base exprssion
%           model and the single repressor model
%   2023-02-02 - Helen Scott
%       * Initial implementation
%       * Copied random_exploration.m and edited to run the variables one-
%           at-a-time

% Load the model catalog
model_catalog;

% Load the parameters
base_parameters;
base_interference_matrix;

parameters = unpackIntMatrix(parameters, i_matrix);

% Get the names of all the parameters because we want to perturb every 
% parameter
parameterNames = keys(parameters);

% Make a dictionary of the description of the V1 value (a string) as the
% key, and the V1 value to use as the value
values = [0, 1];
names = ["Off", "On"];
v1Levels = containers.Map(names, values);

% Set the timespan
tspan = [0 200];

% Set the number of runs to do
nRuns = 5000; 

% Set the output path
results_path = './gRNA-on-off-results/';

% If directory doesn't exist, try to create it
if ~isfolder(results_path)
    disp(['Directory does not exist, attempting to make it: ', ...
        results_path]);
    mkdir(results_path);
end
    
% Loop through just the base expression and the single repressor (using
% their indecies in the model_catalog)
for i=[1,2]
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};

    % Print update on model to screen
    disp(['Currently on: ', modelName]);

    % Clean up name (for making directory/files)
    cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
        ' ', '_'), '/', '-');

    % Set the output path for this model
    outpath = [results_path, cleanModelName, '/'];

    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end

    for k = keys(v1Levels)
        % Get string description of v1 level (for use in plots/files)
        key = k{1};

        % Set the initial
        initial = containers.Map();
        initial('V1') = v1Levels(key);
        initial('V2') = 1;
    
        % Do the perturbation
        results = logNormalPerturbation(modelFun, parameters, ...
            parameterNames, nRuns, tspan, initial, true, true);

        % Extract just what is needed for plotting
%         percentiles = results{:, 1};
%         perturbedValues = results{:, 2};
%         x = results{1, 4}{1};
%         ys = cell(length(results), 1);
%         for k = 1:length(results)
%             ys{k} = results{k, 4}{2};
%         end
    
        % Save results as a .mat file
        save([outpath, 'random-perturb-', cleanModelName, '-gRNA-', key, '.mat'], 'results')
    end
    
end
