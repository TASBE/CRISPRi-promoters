% This file is a list of all models in the order that we'd like to explore and plot them

addpath('../models');

% Each cell row is: {name, function}
models = {
    'Base Expression (No Regulation)'   @No_gRNA_control                    []
    'Single gRNA Repressor'             @Single_gRNA_repression             ["alpha_r_gRNA1"]
    % 'Single gRNA Activator'             @Single_gRNA_activation_traditional ['alpha_r_gRNA1']
    % 'Multiplexed Activation'            @Multiplexed_gRNA_activation        '[alpha_r_gRNA1', 'alpha_r_gRNA2']
    '2 Heterogeneous Target Sites' @Multiplexed_2_gRNA_Repression        ["alpha_r_gRNA1", "alpha_r_gRNA2"]
    '3 Heterogeneous Target Sites' @Multiplexed_3_gRNA_Repression        ["alpha_r_gRNA1", "alpha_r_gRNA2", "alpha_r_gRNA3"]
    '4 Heterogeneous Target Sites' @Multiplexed_4_gRNA_Repression        ["alpha_r_gRNA1", "alpha_r_gRNA2", "alpha_r_gRNA3", "alpha_r_gRNA4"]
    '5 Heterogeneous Target Sites' @Multiplexed_5_gRNA_Repression      ["alpha_r_gRNA1", "alpha_r_gRNA2", "alpha_r_gRNA3", "alpha_r_gRNA4", "alpha_r_gRNA5"]
    '6 Heterogeneous Target Sites' @Multiplexed_6_gRNA_Repression      ["alpha_r_gRNA1", "alpha_r_gRNA2", "alpha_r_gRNA3", "alpha_r_gRNA4", "alpha_r_gRNA5", "alpha_r_gRNA6"]
    '2 Identical Target Sites'    @Multisite_2_gRNA_Repression          ["alpha_r_gRNA1"]
    '3 Identical Target Sites'    @Multisite_3_gRNA_Repression          ["alpha_r_gRNA1"]
    '4 Identical Target Sites'    @Multisite_4_gRNA_Repression          ["alpha_r_gRNA1"]
    '5 Identical Target Sites'    @Multisite_5_gRNA_Repression        ["alpha_r_gRNA1"]
    '6 Identical Target Sites'    @Multisite_6_gRNA_Repression        ["alpha_r_gRNA1"]
    };
n_models = size(models,1);

MODEL_NAME = 1;
MODEL_FUN = 2;
MODEL_PARAMS = 3;
