function rainbowplot(results, variable, variablePos, centerPos, path)
    % FUNCTION NAME:
    %   rainbowplot
    %
    % DESCRIPTION:
    %   Create and save a plot where each row of the result is a different
    %   color line on the plot. Colors are determined by the value of the
    %   perturbed/scanned variable.
    %
    % INPUT:
    %   results - (cell) The output from logNormalPerturbation (Column 1 = 
    %       value perturbed value used in the simulation, column 2 = 
    %       varaibles used for the differential equation, column 3 = 
    %       results from the differential equation (sol))
    %   variable - ()
    %   variablePos - ()
    %   path - (char) output path, will have 'plots/' added to it, default
    %       = './'
    %
    % OUTPUT:
    %   .eps, .fig, .png figure saved in the 'plots/' subfolder of outpath
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   Path is assumed to end in '/' if supplied

    %% Make outout path path
    % Default path is current path
    if nargin < 5
        path = './plots/';  % Note: frontslash works for both Windows and Mac/Unix
    else
        path = [path, 'plots/'];
    end 
    
    % If directory doesn't exist, try to create it
    if ~isfolder(path)
        sanitized_path = strrep(path, '/', '&#47;');
        sanitized_path = strrep(sanitized_path, '\', '&#92;');
        sanitized_path = strrep(sanitized_path, ':', '&#58;');
        fprintf('Directory does not exist, attempting to create it: %s',sanitized_path);
        mkdir(path);
    end
    
    %% Plot
    % Generate colors for all runs
    jetcustom = jet(length(results));

    % Plot
    figure('Visible','off', 'PaperUnits','inches', 'PaperPosition', [0 0 2.75 1.43]);
    for i = 1:length(results)
        tint = results{i, 2};
        yint = results{i, 3};
        if i == centerPos
            plot(tint, yint, 'Color', 'black', 'LineWidth', 2);
        else
            plot(tint, yint, 'Color',  jetcustom(i, :));
        end
        hold on;
    end
    colormap(jet(length(results)));
    
    if length(variable) == 1
        % Color bar for actual variable values
        caxis([results{1,1} results{end,1}]);
        cb = colorbar();
        set(gca,'ColorScale','log')
        % Load the latex version of the parameter and use in color bar
        % label
        parameter_symbols;
        ylabel(cb, sprintf('Value of %s', symbols(variable{1})));
        % title(sprintf('Perturbation of %s', variable{1}));
    elseif length(variable) == 2
        % Color bar for percentiles
        cb = colorbar;
        caxis([results{1, 1} results{end, 1}]);
        ylabel(cb,'Percentile of Perturbed Varriable(s)');
        title(sprintf('Perturbation of %s and %s together', variable(1), variable(2))); 
    else
        cb = colorbar;
        caxis([results{1, 1} results{end, 1}]);
        ylabel(cb,'Percentile of Perturbed Varriable(s)');
        title(sprintf('Joint Perturbation for %s', variable)); % FIXME: Handle multiple names better
    end
    grid on;
    xlabel('Time (Hours)');
    ylabel('Fold Repression');
    set(gca, 'YScale', 'log', 'YLim', [10^0 10^3]);
    fontsize(gca, 10, "points")

    %%  Save figure
    % Output
    saveas(gcf, [path 'scan-', variable{1}, '-fold-repression.png'], 'png');
end
