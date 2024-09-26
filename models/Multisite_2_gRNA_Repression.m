function [time_interval, species_names, y_out, y] = Multisite_2_gRNA_Repression(time_span, parameters, initial, step)
% time_span is the hours values [start, stop]
% parameters is a Map of names to numbers (e.g., rate constants, decay rates, Hill coefficients)
% initial is a Map of variable names to initial values
% step is the number of hours between samples in output; defaults to 1
% Returns vector of time, matrix of output levels at those time points, matrix of all species
    if nargin < 4, step = 1; end
    
    % Define names for input/output variable indexes
    GFP = 1;
	V1 = 2;
	V2 = 3;

    % Set initial values
    y0=zeros(1,6);
    y0(V1) = initial('V1');
	y0(V2) = initial('V2');

    % Set the species names (in the same order as in the ODE)
    species_names = ["GFP", "V1", "V2", "dCas9", "dCas9_gRNA1", "gRNA1"];
    
    % Run ODE
    solution = ode15s(@(t,x) diff_eq(t, x, parameters, species_names), time_span, y0);
    
    % Evaluate species levels at given times
    time_interval = time_span(1):step:time_span(end);
    y = deval(solution, time_interval);
    y_out = y([GFP],:);
end

% ODE differential function
function dx=diff_eq(t, x, parameters, species_names)
    % Unpack parameters from parameter map (and the i_matrix)
    Cas_gRNA_binding = parameters('Cas_gRNA_binding');
	K_R = parameters('K_R');
	alpha_p_GFP = parameters('alpha_p_GFP');
	alpha_p_dCas9 = parameters('alpha_p_dCas9');
	alpha_r_gRNA1 = parameters('alpha_r_gRNA1');
	delta_g = parameters('delta_g');
	lambda = parameters('lambda');
	n = parameters('n');
    

    % Unpack individual species from x
    x = max(1e-12,real(x)); % Truncate values just above zero
    sp = packSpeciesStruct(species_names, x);

    % Compute derivative for each species
    d_V2 = - lambda*sp.V2;
	d_gRNA1 =  alpha_r_gRNA1*sp.V1 - Cas_gRNA_binding*sp.gRNA1*sp.dCas9 - delta_g*sp.gRNA1;
	d_dCas9_gRNA1 =  Cas_gRNA_binding*sp.gRNA1*sp.dCas9 - lambda*sp.dCas9_gRNA1;
	d_dCas9 =  alpha_p_dCas9*sp.V2 - Cas_gRNA_binding*sp.gRNA1*sp.dCas9 - lambda*sp.dCas9;
	d_GFP =  alpha_p_GFP*(K_R^n)/(K_R^n + (sp.dCas9_gRNA1)^n)*(K_R^n)/(K_R^n + (sp.dCas9_gRNA1)^n)*sp.V2 - lambda*sp.GFP;
	d_V1 = - lambda*sp.V1;

    % Pack derivatives for return, ensuring none are complex or go below zero
    dx = max(-x,real([d_GFP, d_V1, d_V2, d_dCas9, d_dCas9_gRNA1, d_gRNA1])');
end
