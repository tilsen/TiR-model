function [] = fbmod_gen_figures()

dbstop if error; close all;
h = fbmod_helpers();

%in order of appearance in manuscript:
fig_functions = {
    'fig_intro';
    'fig_gestural_activation';
    'fig_feedback_concept';
    'fig_timers_actions_autonomy';
    'fig_timers_sens_inter_intra';
    'fig_timers_oscillators';
    'fig_threebody_predictions';
    'fig_CVC_variations';
    'fig_external_influences';
    'fig_competitive_selection';
    'fig_multilevel_competitive_selection';
    'fig_speechrate_variation';
    'fig_projections';
    'fig_prosodic_phrasing'};


for i=1:length(fig_functions)
    eval([fig_functions{i} ';']);
    imfile = dir([h.figures_dir fig_functions{i} '.tif']);
    copyfile([imfile.folder filesep imfile.name],[h.figures_dir sprintf('figure_%02i.tif',i)]);
end


end