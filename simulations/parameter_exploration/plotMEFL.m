function plotMEFL(modelInfo, parameters, parametersToScan, outputFolder)
    %
    % Get the model name and function handle from the model info
    modelName = modelInfo{1};
    cleanModelName = regexprep(modelName, ' ', '_');
    % Get the location where the files are saved
    resultsPath = [outputFolder, '/', cleanModelName, '/'];

    % Define the plots path as the output (normally done in the rainbowplot
    % function)
    outpath = [resultsPath, 'plots/'];
    
    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end

    % Loop through all the parameters
    for parameterName = keys(parametersToScan)
        % Get the position of that parameter in the parameters map
        % Needed because the results just save the values, not the full map
        parameterPos = find(strcmp(keys(parameters), parameterName{1}));
    
        % Get a clean version of the parameter name (e.g. replace underscores
        % with spaces)
        % Keep it in a cell array like parameterName is
        cleanParameterName = strrep(parameterName, '_', ' ');
    
        % TODO: Make a loading function so I don't repeat the code twice
        % Load the gRNA off results
        offFileName = [resultsPath, 'scan-', parameterName{1}, '-gRNA-Off.mat']; % TODO: Don't hardcode the Off
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
        onFileName = [resultsPath, 'scan-', parameterName{1}, '-gRNA-On.mat']; % TODO: Don't hardcode the On
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
    
        % Find the position in the gRNA off results where the base
        % parameter value is used
        offCenterPos = find([gRNA_off_results{:,1}] == parameters(parameterName{1}));
    
        % Find the position in the gRNA off results where the base
        % parameter value is used
        onCenterPos = find([gRNA_on_results{:,1}] == parameters(parameterName{1}));
    
        % Create a list of labels to use
        % Comma separated list of labels to use
        % Only have labels for the base values (thick black lines), for 
        % everything else specify the corresponding label as an empty character
        % vector so that the legend skips that line
        labels = repmat({''}, length(gRNA_off_results) + length(gRNA_on_results), 1);
        for idx = 1:length(gRNA_off_results)
            if idx == offCenterPos
                labels{idx} = "gRNA Off";
            end
        end
        for idx = 1:length(gRNA_on_results)
            if idx == offCenterPos
                labels{length(gRNA_off_results) + idx} = "gRNA On";
            end
        end
    
        
        % TODO: Use a rainbowplot function instead of copying the code and
        % using it twice here
        % Generate colors for all runs
        jetcustom = jet(length(gRNA_on_results)); % Assuming that the on and the off results are the same length
        
        % Make the figure
        figure;
        set(gcf, 'position', [1, 1, 6, 6]);
        hold on;
    
        % Plot the gRNA off rainbow
        for i = 1:length(gRNA_off_results)
            tint = gRNA_off_results{i, 3}{1};
            yint = gRNA_off_results{i, 3}{3};
            if i == offCenterPos
                plot(tint, yint, 'Color', 'black', 'LineWidth', 2, 'LineStyle','--');
            else
                plot(tint, yint, 'Color',  jetcustom(i, :), 'LineStyle','--');
            end
            hold on;
        end
    
        % Plot the gRNA on rainbow
        for i = 1:length(gRNA_on_results)
            tint = gRNA_on_results{i, 3}{1};
            yint = gRNA_on_results{i, 3}{3};
            if i == offCenterPos
                plot(tint, yint, 'Color', 'black', 'LineWidth', 2);
            else
                plot(tint, yint, 'Color',  jetcustom(i, :));
            end
            hold on;
        end
    
        colormap(jet(length(gRNA_on_results)));
        % Color bar for actual variable values
        caxis([gRNA_on_results{1,1} gRNA_on_results{end,1}]);
        cb = colorbar();
        set(gca,'ColorScale','log')
        ylabel(cb, sprintf('Perturbed Value of %s', cleanParameterName{1}));
        title(sprintf('Perturbation of %s', cleanParameterName{1}));
        
        % Style
        grid on;
        xlabel('Time (Hours)');
        ylabel('GFP (MEFL)');
        set(gca, 'YScale', 'log');
        % Add dashed vs solid legend
        legend(labels, 'Location', 'northwest');
    
        % Output
        saveas(gcf, [outpath 'scan-', parameterName{1}, '-mefl.png'], 'png');
    
    end
end