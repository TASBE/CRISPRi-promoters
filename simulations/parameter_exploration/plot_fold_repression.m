% plot_fold_repression
% Last updated: 2023-02-13, by Helen Scott
%
% REVISION HISTORY:
%   2023-02-13 - Helen Scott
%       * Renamed to "plot_fold_repression"
%       * Retitle the plots to be more clear
%   2023-02-10 - Helen Scott
%       * Changed to off/on ratio plotted on a log axis so that the number
%           can be read directly as "Nx repression"
%   2023-02-09 - Helen Scott
%       * Renamed to gRNA_on_off_plots
%       * Change to load the ratio results
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
nBins = [200 50];

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

    % Load the on-off-ratio
    load([modelPath, 'fold-repression.mat'])

    % Convert the lines to points
    xPoints = [];
    yPoints = [];
    for j = 1:length(fold_repression)
        xPoints = [xPoints, x/24];
        yPoints = [yPoints, fold_repression{j}];
    end

    % Plot
    figure('visible', 'off'); % Make but don't show the figure
    hist3([xPoints', yPoints'], nBins, 'CdataMode','auto', 'EdgeColor', 'none');
    xlabel('Time (Days)');
    set(gca, 'YScale', 'log', 'YLim', [1e0 1e1]) % These limits may have to be changed for other runs
    ylabel('Fold Repression');
    cb = colorbar;
    caxis(cLimits);
    ylabel(cb,'Number of Trajectories Per Cell');
    view(2);
    title(sprintf('Fold Repression for Random Perturbations of %s', modelName));

    saveas(gcf, [outpath, 'fold-repression.png'], 'png');

end