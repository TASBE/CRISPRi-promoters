function [n_models_of_interest, models_of_interest, MODEL_NAME, MODEL_FUN] ...
    = collect_models(models_of_interest_names)
    % FUNCTION NAME:
    %   logNormalPerturbation
    %
    % DESCRIPTION:
    %   Generate a series of values for one or more variables perturbed on
    %   a log normal scale. The variables are then passes to a specified
    %   model differential equation.
    %
    % INPUT:
    %   fcnHandle - (function handle) Function handle for differential 
    %       equation function to be used
    %   vars - (map) Map object listing the parameter names and values as 
    %       determined from parameter fitting and literature
    %   perturbedVars - (cell array) Name(s) of variables to be perturbed
    %   nRuns - (double) Number of perturbations to run
    %   tspan - (double) Time span to model ([start, end])
    %   initial - (map) Map object listing initial value names and values
    %   random - (logical) Is the perturbation to be random on the log 
    %       normal scale?
    %   setSeed - (logical) Do you want to set the random number generator
    %       to the default seed (0)? True for yes
    %
    % OUTPUT:
    %   results - (cell) Column 1 = percentile of the perturbed value in
    %       the log normal distribution (NaN for random distributions), 
    %       column 2 = variables used for the differential equation, column
    %       3 = results from the model function
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   Currently asusming std devs for all parameters are 2
    %
    % REVISION HISTORY:
    %   11/14/2021 - Helen Scott
    %       * Change to use with new models
    %           * Got rid of boolean vector (whichVars) because it won't 
    %               work with the parameter map
    %   04/08/2021 - Helen Scott
    %       * Fix for loop (Was recreating varsToUse every loop)
    %

    model_catalog;
    
    % Count how many models to collect
    n_models_of_interest = length(models_of_interest_names);

    % Allot a cell array for the subset of the model catalog I wil want
    models_of_interest = cell(n_models_of_interest, 2);

    % Set a starting index to fill the new model catalog
    model_idx = 1;

    % Look though the full model catalog
    for i = 1:n_models
        % If the model name is in the models_of_interest_names list
        if any(contains(models_of_interest_names, models(2, 1)))
            % Add that row to the models_of_interest catalog
            % For now just taking the name and the function handle, since I
            % haven't been using the variable names at all
            models_of_interest(model_idx, :) = models(i, 1:2);
            % Increase the model index
            model_idx = model_idx + 1;
        end
    end


end