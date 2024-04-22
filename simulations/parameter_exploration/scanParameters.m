function scanParameters(modelInfo, v1Levels, tspan, parameters, parametersToScan, outputFolder)
    % FUNCTION NAME:
    %   scanParameters
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
    
    % Get the model name and function handle from the model info
    modelName = modelInfo{1};
    cleanModelName = regexprep(modelName, ' ', '_');
    modelFun = modelInfo{2};

    % Make an output directory
    outpath = [outputFolder, '/', cleanModelName, '/'];
    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end

    %% Scan through each parameter
    % For each parameter that I want to scan
    for parameterNameCell = keys(parametersToScan)
    
        % Create evenly spaced values for the decades in that range (10 per decade)
        parameterName = parameterNameCell{1}; % Yes, it is always supposed to be 1
        baseParameterValue = parameters(parameterName);
        scanDecades = parametersToScan(parameterName);
        % Log10 value for the low end of the parameter scan
        logLow = log10(baseParameterValue * 10^(-1 * scanDecades(1)));
        % Log10 value for the high end of the parameter scan
        logHigh = log10(baseParameterValue * 10^scanDecades(2));
        % Set the number of values I want (10 per decade plus 1)
        nScans = 10 * sum(scanDecades) + 1;
        scanParameterValues = logspace(logLow, logHigh, nScans);
    
        % Run all the values for this parameter with both the gRNA on and off
        for gRnaKeys = keys(v1Levels)
            % Get string description of v1 level (for use in plots/files)
            key = gRnaKeys{1}; % Yes, it is always supposed to be one
            
            % Display the parameter name and gRNA state so I can keep track
            disp(['Currently on: ', parameterName, '; gRNA ', key])
        
            % Set the initial
            initial = containers.Map();
            initial('V1') = v1Levels(key);
            initial('V2') = 3;
    
            % Make a cell array to hold the results for this paramter and gRNA
            % combination
            results = cell(length(scanParameterValues), 3);
        
            % For every value of that parameter
            for scanIdx = 1:nScans
                % Modulate the parameters
                parametersToUse = containers.Map(parameters.keys, parameters.values);
                parametersToUse(parameterName) = scanParameterValues(scanIdx);
        
                % Run the simulation
                [x, sp, y_out, y] = modelFun(tspan, parametersToUse, initial);
                
                % Add the results to the results variable
                results{scanIdx, 1} = scanParameterValues(scanIdx);
                % Cannot add the parametersToUse as a map, because it was
                % updating to old values to the newest with every loop (i.e. it
                % was not actually saving the correct values for each row, only
                % the last one)
                results{scanIdx, 2} = values(parametersToUse);
                results{scanIdx, 3} = {x + parameters('initial_delay'), sp, y_out, y}; % Shift the x value by the initial lag
    
                % Clear parameters to use?
                clear("parametersToUse")
    
            end
    
            % Save the results for this parameter/gRNA state combination
            save([outpath, 'scan-', parameterName, '-gRNA-', key, '.mat'], 'results')
    
        end
    
    end
end