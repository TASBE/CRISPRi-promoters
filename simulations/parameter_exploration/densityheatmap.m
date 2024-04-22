function densityheatmap(x, ys, nBins, cLimits, name, path)
    % FUNCTION NAME:
    %   densityheatmap
    %
    % DESCRIPTION:
    %   Create ans save a plot of the 
    %
    % INPUT:
    %   results - (cell) The output from logNormalPerturbation (Column 1 = 
    %       percentile of the perturbed value in the log normal 
    %       distribution (NaN for random distributions), column 2 = 
    %       varaibles used for the differential equation, column 3 = 
    %       results from the differential equation (sol))
    %   nBins - (double) The number of bins in each dimension of the 
    %       histogram
    %   cLimits - (double) The colormap limits: a two-element vector of the
    %       form [cmin cmax]
    %   name - (char) Name of the model to be used in plot title
    %   path - (char) output path, will have 'plots/' added to it, default
    %       = './'
    %
    % OUTPUT:
    %   .eps, .fig, .png figure saved in the 'plots/' subfolder of outpath
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   Path is assumed to end in '/' if supplied
    %
    % REVISION HISTORY:
    %   12/06/2021 - Helen Scott
    %       * Changes to work with new results structure
    %   04/08/2021 - hscott
    %       * Initial implement
    %

    %% Add points to arrays
    xPoints = [];
    yPoints = [];
    for i = 1:length(ys)
        xPoints = [xPoints, x/24];
        yPoints = [yPoints, arrayfun(@log10Skip0s, ys{i})];
    end

    %% Plot
    figure('visible', 'off', 'PaperUnits','inches', 'PaperPosition', [0 0 3 1.5]); % Make but don't show the figure
    hist3([xPoints', yPoints'], nBins, 'CdataMode','auto', 'EdgeColor', 'none');
    xlabel('Time (Days)');
    ylim([4.6 7.5])
    ylabel('Log10 [GFP]');
    colorbar('off');
    caxis(cLimits);
%     ylabel(cb,'Number of Trajectories Per Cell');
    view(2);
    title(name);

    % Set the background color to the lowest color in the colormap (to get
    % rid of all the extra white space from setting the ylims)
    colormapCopy = colormap;
    set(gca, 'Color', colormapCopy(1, :))
    
    %%  Save figure
    set(gcf, 'InvertHardcopy', 'off')
    saveas(gcf, [path, 'random-perturb-', name, '.png'], 'png');
end

% Define a function for taking the log10 of something but never return -inf
function logValue = log10Skip0s(linearValue)
    if linearValue == 0
        logValue = nan;
    else
        logValue = log10(linearValue);
    end
end