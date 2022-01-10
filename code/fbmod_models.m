function [M] = fbmod_models(varargin)

%% ---------CVC models
colors = lines(3);
colors = colors([1 2 1 3 3],:);
colors([3 5],:) = pastelize(colors([3 5],:),0.5);

m = fbmod_init_model({'C_1' 'V' 'R_1' 'c_2' 'r_2'},'dur',0.800,'color',colors');

M{1} = max_internal(m); %maximally internal
M{1}.model_descr = {'maximally internal'};

M{end+1} = max_sensory(m); %maximally sensory
M{end}.model_descr = {'maximally sensory'};

M{end+1} = max_oscillators(m); %maximally oscillators
M{end}.model_descr = {'oscillator triggered'};

M{end+1} = hybrid_sensory(m); %
M{end}.model_descr = {'hybrid sensory'};

M{end+1} = hybrid_internal(m); %
M{end}.model_descr = {'hybrid internal'};

M{end+1} = alternate_internal(m); %
M{end}.model_descr = {'alternate internal'};

M{end+1} = hybrid_model(m); %
M{end}.model_descr = {'hybrid model'};

%% ---------CV model
colors = lines(2);
colors = colors([1 2 1],:);
colors(3,:) = pastelize(colors(3,:),0.5);
m = fbmod_init_model({'C' 'V' 'R'},'dur',0.700,'color',colors','dt',0.0001);

M{end+1} = CV_model(m);
M{end}.model_descr = {'CV model'};

%% ---------- gg models
colors = lines(2);
m = fbmod_init_model({'g_1' 'g_2'},'dur',0.500,'color',colors');

M{end+1} = autovsnonauto(m); %
M{end}.model_descr = {'auto vs nonauto'};

m = fbmod_init_model({'g_1' 'g_2'},'dur',0.600,'color',colors');
M{end+1} = externalvsinternal(m); %
M{end}.model_descr = {'external vs internal'};

M{end+1} = intervsintra(m); %
M{end}.model_descr = {'inter vs intra'};

%% ggg models

m = fbmod_init_model({'g_1','g_2','g_3'},'dur',0.500,'freq',5.0,'dt',0.001);
M{end+1} = oscillators_ggg(m); %
M{end}.model_descr = {'oscillators_ggg'};

M{end+1} = external_ggg(m); %
M{end}.model_descr = {'external_ggg'};

M{end+1} = internal_ggg_0(m); %
M{end}.model_descr = {'internal_ggg_0'};

M{end+1} = internal_ggg_1(m); %
M{end}.model_descr = {'internal_ggg_1'};

M{end+1} = internal_ggg_2(m); %
M{end}.model_descr = {'internal_ggg_2'};

M{end+1} = internal_ggg_3(m); %
M{end}.model_descr = {'internal_ggg_3'};

%% select models
if ~isempty(varargin)
    models = varargin{1};
    if ~iscell(models)
        models = {models};
    end
    
    [~,ixocb] = ismember(models,cellfun(@(c)c.model_descr,M));
    
    if any(ixocb==0)
        fprintf('missing model: %s\n', models{ixocb==0}); return;
    end
    
    M = M(ixocb);
    
    if length(M)==1
        M = M{1};
    end
    
end


end

%% oscillators
function [M] = oscillators_ggg(M)

%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3'},...
    {'gate_1','gate_2','gate_3'},0,1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16]+pi;

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3'});
M.phi(osc_ixs,osc_ixs) = 100*[0 1 -1; 1 0 1; -1 1 0];

M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [0 0 0],[1 1 1]);

M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [.075 .075 .075],[-1 -1 -1]);

M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3'},...
    {'gate_1','gate_2','gate_3'},...
    [.075 .195 .075],[-1 -1 -1]);

M = fbmod_remove_systems(M,{'cq','cqg'});

end

%% internal ggg
function [M] = internal_ggg_0(M)
M = fbmod_add_action(M,{'extr_1','extr_2','extr_3'},{'gate_1','gate_2','gate_3'},[0.010 0.050 0.100],1);
M = fbmod_add_action(M,{'etmr_1','etmr_2','etmr_3'},{'gate_1','gate_2','gate_3'},[0.15 0.15 0.15],-1);
end

%% internal ggg
function [M] = internal_ggg_1(M)
M = fbmod_add_action(M,{'etmr_1','etmr_1'},{'gate_2','gate_3'},[0.050 0.100],1);
M = fbmod_add_action(M,{'etmr_1','etmr_2','etmr_3'},{'gate_1','gate_2','gate_3'},[0.15 0.15 0.15],-1);
end

%% internal ggg
function [M] = internal_ggg_2(M)
M = fbmod_add_action(M,{'etmr_1','etmr_2'},{'gate_2','gate_3'},[0.050 0.050],1);
M = fbmod_add_action(M,{'etmr_1','etmr_2','etmr_3'},{'gate_1','gate_2','gate_3'},[0.15 0.15 0.15],-1);
end

%% internal ggg
function [M] = internal_ggg_3(M)
M = fbmod_add_action(M,{'etmr_1','extr_2'},{'gate_3','gate_2'},[0.050 0.050],1);
M = fbmod_add_action(M,{'etmr_1','etmr_2','etmr_3'},{'gate_1','gate_2','gate_3'},[0.15 0.15 0.15],-1);
end

%% external ggg
function [M] = external_ggg(M)

M = fbmod_add_action(M,...
    {'stmr_1','stmr_2'},...
    {'gate_2','gate_3'},...
    [0.050 0.05],[1 1]);

M = fbmod_add_action(M,...
    {'stmr_1','stmr_2','stmr_3'},...
    {'gate_1','gate_2','gate_3'},...
    [0.1 0.1 0.1],[-1 -1 -1]);

end

%% mixed ggg
function [M] = mixed_ggg(M)

M = fbmod_add_action(M,...
    {'etmr_1','etmr_1'},...
    {'gate_2','gate_3'},...
    [0.050 0.075],[1 1]);

M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3'},...
    {'gate_1','gate_2','gate_3'},...
    [0.1 0.1 0.1],[-1 -1 -1]);

end

%% inter vs. intra
function [M]= intervsintra(M)

M = fbmod_add_action(M,'itmr_1','gate_1',0.250,-1);
M = fbmod_add_action(M,'etmr_1','gate_2',0.200,1);
M = fbmod_add_action(M,'etmr_2','gate_2',0.200,-1);

end

%% external vs internal
function [M]= externalvsinternal(M)

M = fbmod_add_action(M,'etmr_1','gate_1',0.250,-1);
M = fbmod_add_action(M,'etmr_1','gate_2',0.250,1);
M = fbmod_add_action(M,'stmr_2','gate_2',0.200,-1);

end

%% auto vs. non-auto
function [M]= autovsnonauto(M)
M = fbmod_add_action(M,'extr_2','gate_1',0.250,-1);
M = fbmod_add_action(M,'extr_1','gate_2',0.200,1);
M = fbmod_add_action(M,'etmr_2','gate_2',0.250,-1);
M = fbmod_remove_systems(M,{'cq','cqg','osc','oscr','oscr2','oscg'});
end

%% hybrid sensory
function [M]= hybrid_sensory(M)

M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3'},...
    {'gate_1','gate_2','gate_3'},...
    [0 0 0],1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16 0 0]+pi;

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3'});
M.phi(osc_ixs,osc_ixs) = 100*...
    [ 0  1 -1;
      1  0  1; 
     -1  1  0];
 
%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

%termination
M = fbmod_add_action(M,...
    {'etmr_1' 'stmr_2' 'etmr_3' 'stmr_4' 'stmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [0.100 0.250 0.100 0.100 0.100],-1);  

%post-vocalic initiation
M = fbmod_add_action(M,...
    {'stmr_2','stmr_4'},...
    {'gate_4','gate_5'},...
    [0.250 0.100],1);  %C1 supression

M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [0 0 0],[1 1 1 1 1]);

M = fbmod_add_action(M,...
    {'extr_2','extr_2','extr_2'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [.2 .2 .2],[-1 -1 -1 -1 -1]);

end

%% hybrid internal
function [M]= hybrid_internal(M)

M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3'},...
    {'gate_1','gate_2','gate_3'},...
    [0 0 0],1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16 0 0]+pi;

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3'});
M.phi(osc_ixs,osc_ixs) = 100*...
    [ 0  1 -1;
      1  0  1; 
     -1  1  0];
 
%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

%termination
M = fbmod_add_action(M,...
    {'etmr_1' 'etmr_2' 'etmr_3' 'etmr_4' 'etmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [0.100 0.250 0.100 0.100 0.100],-1);  

%post-vocalic initiation
M = fbmod_add_action(M,...
    {'etmr_2','etmr_4'},...
    {'gate_4','gate_5'},...
    [0.200 0.075],1);  %C1 supression

M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [0 0 0],[1 1 1 1 1]);

M = fbmod_add_action(M,...
    {'extr_2','extr_2','extr_2'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [.2 .2 .2],[-1 -1 -1 -1 -1]);

end

%% hybrid model
function [M]= hybrid_model(M)

M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3'},...
    {'gate_1','gate_2','gate_3'},...
    [0 0 0],1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16 0 0]+pi;

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3'});
M.phi(osc_ixs,osc_ixs) = 100*...
    [ 0  1 -1;
      1  0  1; 
     -1  1  0];
 
%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3','etmr_4','etmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [0.050 0.250 0.050 0.050 0.050],-1);  %V activation

M = fbmod_add_action(M,...
    {'etmr_4'},...
    {'gate_5'},...
    [0.050],[1]);  %C1 supression

M = fbmod_add_action(M,...
    {'etmr_2','etmr_2'},...
    {'gate_2','gate_4'},...
    [0.250 0.200],[-1 1]);  %C1 supression

M = fbmod_add_action(M,...
    {'stmr_2','stmr_2','stmr_4','stmr_4'},...
    {'gate_2','gate_4','gate_4','gate_5'},...
    [0.250 0.200 0.05 0.05],[-1 1 -1 1]);  %C1 supression

M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [0 0 0],[1 1 1 1 1]);

M = fbmod_add_action(M,...
    {'extr_2','extr_2','extr_2'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [.2 .2 .2],[-1 -1 -1 -1 -1]);

end

%% CV model
function [M]= CV_model(M)

M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3'},...
    {'gate_1','gate_2','gate_3'},...
    [0 0 0],1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16]+pi;

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3'});
M.phi(osc_ixs,osc_ixs) = 50*...
    [ 0  1 -1;
      1  0  1; 
     -1  1  0];
 
%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

gdurs = [0.100 0.300 0.100];

M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3'},...
    {'gate_1','gate_2','gate_3'},...
    gdurs,-1);  %V activation

M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [0.01 0.02 0.03],[1 1 1]);

M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3'},...
    {'oscg_1','oscg_2','oscg_3'},...
    0.075,[-1 -1 -1]);

M = fbmod_remove_systems(M,{'cq','cqg'});

end

%% alternate internal
function [M]= alternate_internal(M)

M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3','osc_4','osc_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [0 0 0 0 0],1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16 [pi/16 -pi/16]]+pi;

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3','osc_4','osc_5'});
M.phi(osc_ixs,osc_ixs) = 100*...
    [ 0  1 -1  0  0;
      1  0  1  0  0; 
     -1  1  0  0  0;
      0  0  0  0 -1;
      0  0  0 -1  0];
 
%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

sup_times = [0.050     0.250    0.050   0.050   0.050];

%gestural supression
M = fbmod_add_action(M,...
    {'itmr_1','etmr_2','itmr_3','itmr_4','itmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},sup_times,-1); 

%
M = fbmod_add_action(M,...
    {'itmr_1','itmr_2','itmr_3','itmr_4','itmr_5'},...
    {'oscg_1','oscg_2','oscg_3','oscg_4','oscg_5'},sup_times/2,-1); 


%internal supression of vowel
M = fbmod_add_action(M,...
    {'etmr_2'},...
    {'gate_2'},...
    [0.250],[-1]);  %C1 supression

M = fbmod_add_action(M,...
    {'etmr_2','etmr_2'},...
    {'oscg_4','oscg_5'},...
    [0.1 0.1],[1 1]);  %C1 supression

M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3'},...
    [0 0 0],[1 1 1 1 1]);

end

%% max_internal
function [M]= max_internal(M)

cdur = 0.075;
vdur = 0.300;
ofac = 2;

%initiation
M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_2','etmr_4'},...
    {'gate_2','gate_3','gate_4','gate_5'},...
    [cdur/ofac cdur/ofac vdur cdur/ofac],1); 

%termination
M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3','etmr_4','etmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [cdur vdur cdur cdur cdur],-1);  

end

%% max_oscillators
function [M]= max_oscillators(M)

cdur = 0.075;
vdur = 0.300;
ofac = 2;

%oscillator triggering
M = fbmod_add_action(M,...
    {'osc_1','osc_2','osc_3','osc_4','osc_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [0 0 0 0 0],1);

%initial oscillator phases
M.X0(ismember(M.S.type,'osc')) = [pi/16 0 -pi/16 [pi/16 -pi/16]+pi/2]+pi;

M.omega(ismember(M.S.type,'osc')) = [5 5 5 8 8]*(2*pi);

%oscillator relative phase coupling
osc_ixs = M.getid({'osc_1','osc_2','osc_3','osc_4','osc_5'});
M.phi(osc_ixs,osc_ixs) = 100*...
    [ 0  1 -1  0  0;
      1  0  1 -1  0; 
     -1  1  0  0  0;
      0 -1  0  0 -1;
      0  0  0 -1  0];
  
%termination
M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3','etmr_4','etmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [cdur vdur cdur cdur cdur],-1);  

%oscilator de-gating
M = fbmod_add_action(M,...
    {'etmr_1','etmr_2','etmr_3','etmr_4','etmr_5'},...
    {'oscg_1','oscg_2','oscg_3','oscg_4','oscg_5'},...
    [cdur vdur cdur cdur cdur],-1);

%remove extrinsic trigger
M.chi(M.getid('extr_1'),:) = 0;
M.T = M.T(~ismember(M.T.src,M.getid('extr_1')),:); 

%gate oscillators
M = fbmod_add_action(M,...
    {'extr_1','extr_1','extr_1','extr_1','extr_1'},...
    {'oscg_1','oscg_2','oscg_3','oscg_4','oscg_5'},...
    [0 0 0 .35 .35],1);

end

%% max_sensory
function [M]= max_sensory(M)

cdur = 0.075;
vdur = 0.300;
ofac = 2;

%initiation
M = fbmod_add_action(M,...
    {'stmr_1','stmr_2','stmr_2','stmr_4'},...
    {'gate_2','gate_3','gate_4','gate_5'},...
    [cdur/ofac cdur/ofac vdur cdur/ofac],1);  

%termination
M = fbmod_add_action(M,...
    {'stmr_1','stmr_2','stmr_3','stmr_4','stmr_5'},...
    {'gate_1','gate_2','gate_3','gate_4','gate_5'},...
    [cdur vdur cdur cdur cdur],-1);  
end

