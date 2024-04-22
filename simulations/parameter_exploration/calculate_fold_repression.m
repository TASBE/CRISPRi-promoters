% calculate_fold_repression
% Load the results from the gRNA_on_off perturbations and take a ratio
% between the two states at each time point for each run
% Save the ratio as the "fold repression"
% Last updated: 2023-02-13, by Helen Scott
%
% REVISION HISTORY:
%   2023-02-13 - Helen Scott
%       * Changed to fold repression
%   2023-02-10 - Helen Scott
%       * Changed to off/on ratio

% Load the model catalog
model_catalog;

% Set the v1Levls key names (has to match what was used in gRNA_on_off)
v1LevelNames = ["Off", "On"];

% Set the results path
resultsPath = './gRNA-on-off-results/';

for i=[1,2]
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};

    % Print update on model to screen
    disp(['Currently on: ', modelName]);

    % Clean up name (for making directory/files)
    cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
        ' ', '_'), '/', '-');

    % Set the output path for this model
    modelPath = [resultsPath, cleanModelName, '/'];

    % Load the off results (Assumes that off is always first in the levels
    % list)
    resultsFile = [modelPath, 'random-perturb-', cleanModelName, '-gRNA-', v1LevelNames{1}, '.mat'];
    load(resultsFile)

    % Rename the results variable
    gRNA_off_results = results;
    clear("results")

    % Load the on results
    resultsFile = [modelPath, 'random-perturb-', cleanModelName, '-gRNA-', v1LevelNames{2}, '.mat'];
    load(resultsFile)

    % Rename the results variable
    gRNA_on_results = results;
    clear("results")

    % Extract the time points
    % Only have to do this for one of the simulations because the time
    % points are the same
    x = gRNA_on_results{1, 4}{1};

    % Take the ratio of the runs
    fold_repression = cell(length(gRNA_on_results), 1); % Could be using either set of results
    for k = 1:length(gRNA_on_results)
        fold_repression{k} = 10.^(log10(gRNA_off_results{k, 4}{2}) - log10(gRNA_on_results{k, 4}{2}));
    end

    % Save the ratios in a .mat file
    save([modelPath, 'fold-repression.mat'], 'x', 'fold_repression')

end