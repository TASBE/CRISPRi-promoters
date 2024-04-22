% gRNA_on_off_plots
% Last updated: 2023-02-09, by Helen Scott
%
% REVISION HISTORY:
%   2023-02-09 - Helen Scott
%       * Renamed to gRNA_on_off_plots
%       * Change paths to match where the results are
%   2023-02-02 - Helen Scott
%       * Copied random_exploration_plots.m and made the necessary changes
%           for it to run on the one-at-a-time analysis
%

% Load the model catalog
model_catalog;

% Load the parameters % FIXME: I don't think I need these
base_parameters;
base_interference_matrix;

parameters = unpackIntMatrix(parameters, i_matrix);

% Set the v1Levls key names (has to match what was used in gRNA_on_off)
v1LevelNames = ["Off", "On"];

% Set the results path
resultsPath = './gRNA-on-off-results/';

% Set color bar limits for all plots
cLimits = [0 500];

% Set number of bins for all plots
nBins = [75 25];

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
    modelPath = [resultsPath, cleanModelName, '/'];

    % Make an output folder for the plots for this model
    outpath = [modelPath, 'plots/'];

    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end

    % Loop thorugh the v1 levels used
    for v1LevelIdx = 1:length(v1LevelNames)
        v1Level = v1LevelNames(v1LevelIdx);
        
        % Generate the name of the results file
        resultsFile = [modelPath, 'random-perturb-', cleanModelName, '-gRNA-', v1Level{1}, '.mat'];

        % Load the results
        load(resultsFile)

        % Pull out the xs and the ys
        x = results{1, 4}{1};
        ys = cell(length(results), 1);
        for k = 1:length(results)
            ys{k} = results{k, 4}{2};
        end

        % Plot
        densityheatmap(x, ys, nBins, cLimits, [modelName, ' (gRNA ', v1Level{1}, ')'], outpath);
    end

end