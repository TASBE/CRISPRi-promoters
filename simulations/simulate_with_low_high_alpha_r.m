% Load the parameters and model info
base_parameters;
model_catalog;
base_interference_matrix;

% Set inital values for the simulations
initial = containers.Map();
initial('V1') = 10;
initial('V2') = 10;

% Is it reasonable to adjust the alpha values up and down by 100x?
tuning = [-2 2];
n_tunings = numel(tuning);

time = [0 100];
y_out = zeros(n_models,n_tunings,n_tunings,101);
y_complete = cell(n_models,n_tunings,n_tunings);
for i=1:n_models
    fprintf('Tuning model %s',models{i,MODEL_NAME});
    count = 0;
    for t1=1:n_tunings
        for t2=1:n_tunings
            tuned_param = containers.Map(parameters.keys, parameters.values);
            if models{i,MODEL_PARAM1}, tuned_param(models{i,MODEL_PARAM1}) = tuned_param(models{i,MODEL_PARAM1})*10^tuning(t1); end;
            if models{i,MODEL_PARAM2}, tuned_param(models{i,MODEL_PARAM2}) = tuned_param(models{i,MODEL_PARAM2})*10^tuning(t2); end;
            try
                [time_interval, y_out(i,t1,t2,:), y_complete{i,t1,t2}] = models{i,MODEL_FUN}(time,tuned_param,i_matrix,initial,1);
                count = count+1;
                if mod(count,10)==0, fprintf('.'); end;
            catch
                fprintf('!');
            end
        end
    end
    fprintf('\n');
end

% Not saving this because it ends up being too big, at 68 MB
%save('parameter_exploration.mat','tuning','time_interval','y_out','y_complete');

line_types = {"--", "-.", ":", "-o", '-'};
color = ['b', 'r'];
for i=1:n_models
    h = figure('PaperPosition',[1 1 6 6]); 
    plot(time_interval, log10(squeeze(y_out(1,1,1,:))), 'DisplayName', 'Constitutive'); hold on;
    for t1=1:n_tunings
        if t1==1, t1_description = 'Low'; else, t1_description = 'High'; end
%             line_spec = line_types{t1};
%             color = [1-(t1/n_tunings), 0, t2/n_tunings];
%             plot(time_interval, log10(squeeze(y_out(i,t1,1,:))), ...
%                 line_spec, ...
%                 'Color', color, ...
%                 'DisplayName', ['\alpha_{r,gRNA1} = ', t1_description]);
        for t2=1:n_tunings
            if t2==1, t2_description = 'Low'; else, t2_description = 'High'; end
            color = [1-(t1/n_tunings), 0, t2/n_tunings];
            % Red = low Param 1; Blue = high Param 1; Green = low Param 2
%             color = [max(0,-tuning(t1)/2), 1-t2/n_tunings, max(0,tuning(t1)/2)];
            line_spec = line_types{(t1*2 - 1) + (t2-1)};
            plot(time_interval, log10(squeeze(y_out(i,t1,t2,:))), ...
                line_spec, ...
                'Color', color, ...
                'DisplayName', ['\alpha_{r,gRNA1} = ', t1_description, '; \alpha_{r,gRNA2} = ', t2_description]);
            hold on;
        end
    end
    xlabel('Hours'); ylabel('Log10 [GFP]');
    xlim([0 48]); % ylim([4 7]);
    title(models{i,MODEL_NAME});
    % Add legend
    legend('Location','SouthOutside');
    % Save figure as 3 different file types (uses TASBE)
    outputfig(h,sprintf('low_vs_high_gRNA_%s',models{i,MODEL_NAME}),'plots/low_vs_high_gRNA');
end
