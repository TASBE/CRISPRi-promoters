function plotFoldRepression(modelInfo, parameters, parametersToScan, outputFolder)
    %    
    % Get the model name and function handle from the model info
    modelName = modelInfo{1};
    cleanModelName = regexprep(modelName, ' ', '_');
    % Get the location where the files are saved
    resultsPath = [outputFolder, '/', cleanModelName, '/'];
    
    % Loop through all the parameters
    for parameterName = keys(parametersToScan)
        % Load the fold repression results
        resultsFileName = [resultsPath, 'scan-', parameterName{1}, '-fold-repression.mat'];
        if isfile(resultsFileName)
            load(resultsFileName);
        else
            warning(['No parameter scan was performed for ', ...
                parameterName{1}, '. Skipping.'])
            continue
        end
    
        % Get the position of that parameter in the parameters map
        % Needed because the results just save the values, not the full map
        parameterPos = find(strcmp(keys(parameters), parameterName{1}));
    
        % Find the position in the fold repression results where the base
        % parameter value is used
        centerPos = find([fold_repression{:,1}] == parameters(parameterName{1}));
    
        % Call the rainbowplot function to make an save a plot
        rainbowplot(fold_repression, parameterName, parameterPos, ...
            centerPos, resultsPath);
    
    end

end