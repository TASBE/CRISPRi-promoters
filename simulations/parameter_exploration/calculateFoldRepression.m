function calculateFoldRepression(modelInfo, parametersToScan, outputFolder)
    % Get the location where the files are saved
    % Get the model name and function handle from the model info
    modelName = modelInfo{1};
    cleanModelName = regexprep(modelName, ' ', '_');
    outpath = [outputFolder, '/', cleanModelName, '/'];

   % Loop through all the parameters
    for parameterName = keys(parametersToScan)
        % Load the gRNA off results
        offFileName = [outpath, 'scan-', parameterName{1}, '-gRNA-Off.mat'];
        if isfile(offFileName)
            load(offFileName);
            % Rename the results variable to signify the gRNA state
            gRNA_off_results = results;
            clear("results")
        else
            warning(['No parameter scan was performed for ', ...
                parameterName{1}, '. Skipping.'])
            continue
        end
    
        % Load the gRNA off results
        onFileName = [outpath, 'scan-', parameterName{1}, '-gRNA-On.mat'];
        if isfile(onFileName)
            load(onFileName);
            % Rename the results variable to signify the gRNA state
            gRNA_on_results = results;
            clear("results")
        else
            warning(['No parameter scan was performed for ', ...
                parameterName{1}, '. Skipping.'])
            continue
        end
    
        % Extract the timepoints of the simulation
        % Only have to do this for one of the simulations because the time
        % points are the same
        x = gRNA_on_results{1, 3}{1};
    
        % Take the ratio of the gRNA states
        fold_repression = cell(length(gRNA_on_results), 2); % Could be using either set of results
        for k = 1:length(fold_repression)
            fold_repression{k, 1} = gRNA_off_results{k, 1}; % On and off should have the same value, but should add a check
            fold_repression{k, 2} = gRNA_off_results{k, 3}{1};
            fold_repression{k, 3} = gRNA_off_results{k, 3}{3} ./ gRNA_on_results{k, 3}{3};
        end
    
        % Save the ratios in a .mat file
        save([outpath, 'scan-', parameterName{1}, '-fold-repression.mat'], ...
            'x', 'fold_repression')
    end
end