%% Set-Up
% Set the output folder
outpath = 'vector_concentration_test_results/';

% If directory doesn't exist, try to create it
if ~isfolder(outpath)
    disp(['Directory does not exist, attempting to make it: ', ...
        outpath]);
    mkdir(outpath);
end

% Load the models
models = containers.Map();
models('multisite_gRNA_repression') = @Multisite_gRNA_Repression;
models('multiplexed_gRNA_epression') = @Multiplexed_gRNA_Repression;

% Load the base parameters
base_parameters;

% Not sure if I need the interference matrix at all
base_interference_matrix;
parameters = unpackIntMatrix(parameters, i_matrix);

% Settings for the model simulations
% Set the timespan
tspan = [0 200];

% List of vector concnetrations of interest
vector_concentrations = [1, 10, 100];

%% Run the models
% Loop through each of the models
for modelName = keys(models)
    % Get the model function handle
    modelFun = models(modelName{1}); 

    % Make a cell array for the results
    results = cell(length(vector_concentrations), 2);

    % Loop through the vector concentration values of interest
    for concIdx = 1:length(vector_concentrations)
        % Get the concentration from the list based on the index
        conc = vector_concentrations(concIdx);

        % Save the vector concentration in the results array
        results{concIdx, 1} = conc;
        
        % Set the initial
        initial = containers.Map();
        initial('V1') = conc;
        initial('V2') = conc;

        % Run the model
        [x, y_out, y] = modelFun(tspan, parameters, initial);

        % Save the results to the holder
        results{concIdx, 2} = {x + parameters('initial_delay'), y_out, y}; % Shift the x value by the initial lag
    end

    % Save the holder for this model
    save([outpath, 'scan-vector-concentrations-', modelName{1}, '.mat'], ...
        'results')

end


