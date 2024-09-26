parameters = containers.Map();

% This file contains all of the base parameter values
% These parameters were taken from the NIH-CRISPR work, ideally I would
% like to fit the parameters to a more relevant dataset
% This should be the starting point for any simulation run

% Transcription/translation parameters
parameters('alpha_r_gRNA2') = 10^(3.3090); % Molecules/hour; Taken from NIH-CRISPR
parameters('alpha_r_gRNA1') = parameters('alpha_r_gRNA2');
parameters('alpha_r_gRNA3') = parameters('alpha_r_gRNA2');
parameters('alpha_r_gRNA4') = parameters('alpha_r_gRNA2');
parameters('alpha_r_gRNA5') = parameters('alpha_r_gRNA2');
parameters('alpha_r_gRNA6') = parameters('alpha_r_gRNA2');
parameters('alpha_p_GFP') = 10^(5.5793); % Molecules/hour; Fit to the 2023-03-14 timecourse
parameters('alpha_p_dCas9') = 10^(3.0415); % Molecules/hour; Taken from NIH-CRISPR

% Degradation/dilution parameters
parameters('delta_g') = 10^(0.0003); % Fraction/hour; Taken from NIH-CRISPR; assuming identical for all gRNAs
parameters('lambda') = 10^(-1.6225); % Fraction/hour; Fit to the 2023-03-14 timecourse

% dCas9 activation/repression mechanism parameters
parameters('Cas_gRNA_binding') = 10^(-4.2577);  % Taken from NIH-CRISPR
parameters('K_A') = 2.34 * 10^6; % Used in NIH-CRISPR, from Calin Belta paper; Cannot be readily modulated
parameters('K_R') = 10^2.9; % Hypothesized; Adjusted to give expected fold repression value
parameters('n') = 0.92; % Used in NIH-CRISPR, from Calin Belta paper; Cannot be readily modulated

% Initial delay
parameters('initial_delay') = 18.0149; % Hours; Fit to the 2023-03-14 timecourse
