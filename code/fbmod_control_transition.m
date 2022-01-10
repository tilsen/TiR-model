function [] = fbmod_control_transition()

dbstop if error;
h = fbmod_helpers();

%% determine alphas

lambdas = linspace(0,1,21);

b_etmr = 2;
b_stmr = 1;

   
alpha_etmr = 1 ./ (1+b_etmr*lambdas);
alpha_stmr = 1 ./ (1+b_stmr*lambdas);


%% specify and run models
m = fbmod_init_model({'V_1','c_1'},'dur',0.700,'dt',0.001);
m = fbmod_add_action(m,...
    {'etmr_1','etmr_1','stmr_1','etmr_2'},...
    {'gate_1','gate_2','gate_2','gate_2'},...
    [0.150 0.1 0.1 0.1],...
    [-1 1 1 -1]);


etmr_ix = ismember(m.S.type,'etmr');
stmr_ix = ismember(m.S.type,'stmr');

for i=1:length(lambdas)
  
    mx = m;
    mx.alpha(:,etmr_ix) = mx.alpha(:,etmr_ix) ./ (1+b_etmr*lambdas(i));
    mx.alpha(:,stmr_ix) = mx.alpha(:,stmr_ix) ./ (1+b_stmr*lambdas(i));
    
    M{i} = mx;
    
    R(i,1).model = sprintf('lambdavar_%02i',i);
    R(i).lambda = lambdas(i);
    R(i).alpha_etmr = alpha_etmr(i);
    R(i).alpha_stmr = alpha_stmr(i);
end

M = fbmod_run_models(M);

%% extract gestural timing

M = fbmod_prep_schema(M);
gest_ix = ismember(M{1}.S.type,'gest');
for i=1:length(M)
    [~,~,R(i).t_on,R(i).t_off] = fbmod_gestural_activation(M{i}.X(:,gest_ix),M{i}.t);
end
R = struct2table(R);
R.t_on = cell2mat(R.t_on);
R.t_off = cell2mat(R.t_off);
R.deltas = diff(R.t_on,[],2);


for i=1:length(M)
       
    %
    T = M{i}.T;    
    ix_etmr_force = find(T.F{h.getix(T,'etmr_1')}(:,2)>0,1,'first');
    ix_stmr_force = find(T.F{h.getix(T,'stmr_1')}(:,1)>0,1,'first');
    
    if isempty(ix_stmr_force)
       R.stmr_first(i) = -1;
    elseif ix_stmr_force<ix_etmr_force
       R.stmr_first(i) = 1;
    elseif ix_stmr_force>ix_etmr_force
       R.stmr_first(i) = -1;
    end
    
end

ixs = [1 find(R.stmr_first>0,1,'first') length(M)];
M = M(ixs);

%%
save([h.sim_dir mfilename '.mat'],'R','M');


end

