function [] = fbmod_prosodic_phrasing()

dbstop if error; close all;
addpath('..\model');
h = fbmod_helpers;

fs = [30 24 16 12];

%% ------- prosodic phrasing example
%Al and Bo or Cam were there (each row is a pword)

utterance = {...
    {{'V_1' 'c_1'}}...
    {{'V_2' 'c_{21}' 'c_{22}'} {'C_3' 'V_3'}}...
    {{'V_4' 'c_4'} {'C_5' 'V_5' 'c_5'}}...
    {{'C_6' 'V_6'} {'C_7' 'V_7' 'c_7'}}};

%%

[M,U,P] = fbmod_gen_utterance_model('standard',utterance,'dur',3.5);
M.color = viridis(M.Ng)';
M = fbmod_run_models(M);

%%
schema_params = {...
    'radius',[1/3 1/8 1/5],...
    'fontsize',fs([2 4 3]),...
    'labeltype',{'full' 'minimal' 'full'}};

M = fbmod_prep_schema(M,schema_params{:});

M.PHI_XX = []; M.F_XX = [];

save([h.sim_dir mfilename '.mat'],'M','U','P');


end
