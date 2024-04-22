% arrange_figures
% This is a script that collects the parameter scan plots for a few models
% and arranges them in a grid and then saves that multifigure
% NOTE: The arranging works, but the resolution is super low, so the figure
% is not usable.

% Set the models of interest
% The names used here must match the names in the model catalog exactly
modelsOfInterest = {'Single_gRNA_Repressor', '4_Heterogeneous_Target_Sites', '4_Identical_Target_Sites'};

% Set the parameters to display
parametersToShow = {'Cas_gRNA_binding', 'alpha_p_GFP', 'alpha_p_dCas9', 'delta_g', 'lambda', 'K_R'};

% Read the images into a cell array
imgArray = cell(length(parametersToShow), length(modelsOfInterest));
for modelIdx = 1:length(modelsOfInterest)
    for paramIdx = 1:length(parametersToShow)
        img = imread(['results/', modelsOfInterest{modelIdx}, '/plots/scan-', parametersToShow{paramIdx}, '-mefl.png']);
        imgArray{paramIdx, modelIdx} = img;
    end
end

% Rearragne thecell array to match the montage order (column-major)
montageOrder = reshape(imgArray.', 1, []);

% Resize the images to a specific 


% View the images as a montage
montageObj = montage(montageOrder, 'Size', [6, 3], 'ThumbnailSize', [1427, 2671]);

% Save the montage
imwrite(getframe(gca).cdata, 'results/single_4hetero_4identical.png')