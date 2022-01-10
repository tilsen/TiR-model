function [] = fbmod_speechrate_boundary_variation()

dbstop if error; %close all;
h = fbmod_helpers;

%: model of speech rate variation effects on three word sequence
%: third word has heightened effects


%% ------- utterance model

sn = @(str,n)sprintf('%s_%i',str,n);

sysnames = {...
    {sn('C',1),sn('V',1),sn('R',1),sn('c',1),sn('r',1)};
    {sn('C',2),sn('V',2),sn('R',2)};
    {sn('C',3),sn('V',3),sn('R',3),sn('c',3),sn('r',3)}};

utterance = {sysnames}; %one pwrd

%[M0,S] = fbmod_gen_pword_model('standard',sysnames,'dur',4.00,'color',viridis(numel([sysnames{:}]))'); %#ok<*ASGLU>
[M0,S] = fbmod_gen_utterance_model('standard',utterance,'dur',4.00,'color',viridis(numel([sysnames{:}]))'); %#ok<*ASGLU>

%M0 = fbmod_run_models(M0); fbmod_quickplot(M0);

Nm = 11;


% lambda->omega mapping:
frq_range = 8*(2*pi);
lambda_resp = 5;
base_frq = 6*(2*pi);
frq_eff = @(lambda)base_frq - (frq_range/2)*tanh(lambda_resp*(lambda-1/2));

lambdas = linspace(0,1,Nm);
omegas = frq_eff(lambdas);

b_etmr = 2;
b_stmr = 1;

b_etmr_last = 2.2;
b_stmr_last = 1.2;

gest_ix = ismember(M0.S.type,'gest');
etmr_ix = ismember(M0.S.type,'etmr');
stmr_ix = ismember(M0.S.type,'stmr');
osc_ix = ismember(M0.S.type,'osc');
    
alpha_etmr = 1 ./ (1+b_etmr*lambdas);
alpha_stmr = 1 ./ (1+b_stmr*lambdas);

last_word_ix = ismember(M0.S.type_id,S.gest_typeid{end});

M = {};
for i=1:Nm
  
    m = M0;
    
    %adjust etmr & stmr alphas
    m.alpha(:,etmr_ix & ~last_word_ix) = m.alpha(:,etmr_ix & ~last_word_ix) ./ (1+b_etmr*lambdas(i));
    m.alpha(:,stmr_ix & ~last_word_ix) = m.alpha(:,stmr_ix & ~last_word_ix) ./ (1+b_stmr*lambdas(i));
    
    m.alpha(:,etmr_ix & last_word_ix) = m.alpha(:,etmr_ix & last_word_ix) ./ (1+b_etmr_last*lambdas(i));
    m.alpha(:,stmr_ix & last_word_ix) = m.alpha(:,stmr_ix & last_word_ix) ./ (1+b_stmr_last*lambdas(i));    
    
    %adjust frequencies
    m.omega(osc_ix) = omegas(i);
    
    M{i} = m;
    
    R(i,1).model = sprintf('ratevar_%02i',i);
    R(i).lambda = lambdas(i);
    R(i).omega = omegas(i);
    R(i).alpha_etmr = 1/(1+b_etmr*lambdas(i));
    R(i).alpha_stmr = 1/(1+b_stmr*lambdas(i));
    
    R(i).alpha_etmr_last = 1/(1+b_etmr_last*lambdas(i));
    R(i).alpha_stmr_last = 1/(1+b_stmr_last*lambdas(i));    
end

batchsize = 4;
batch_ixs = {};
for i=1:batchsize:Nm
    batch_ixs{end+1} = [i:min(Nm,i+batchsize-1)]; %#ok<AGROW>
end

for i=1:length(batch_ixs)
    
    tic
    Mbatch = fbmod_run_models(M(batch_ixs{i}),'parallel',true);
    toc
    
    for j=1:length(Mbatch)
        mix = batch_ixs{i}(j);
        [~,~,R(mix).t_on,R(mix).t_off] = fbmod_gestural_activation(Mbatch{j}.X(:,gest_ix),Mbatch{j}.t);
    end
    
    if i==1
        M0 = Mbatch{1};
    end

end
       
%%
R = struct2table(R);
save([h.sim_dir mfilename '.mat'],'R','M0');

end
