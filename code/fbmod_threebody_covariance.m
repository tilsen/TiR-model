function [] = fbmod_threebody_covariance()

dbstop if error; %close all;
addpath('..\model');
h = fbmod_helpers;

%: variance-covariance predictions of various models of three-gesture
%systems, with parametrically varied local and global noise in
%oscillator frequency and growth rates

models = {
    'oscillators_ggg'
    'external_ggg' 
    'internal_ggg_0'
    'internal_ggg_1'
    'internal_ggg_2'
    'internal_ggg_3'}';

noiseparset = [1 2 2 2 2 2];

N_models = length(models);
M = fbmod_models(models);

N_reps = 400;
batch_size = 400;
N_batches = N_reps/batch_size;

%% 

max_omega_range = 2*pi*5;
max_alpha_sigma = 0.2;

ix_osc = ismember(M{1}.S.type,'osc')';
ix_etmr = ismember(M{1}.S.type,'etmr')';
ix_stmr = ismember(M{1}.S.type,'stmr')';
ix_gest = ismember(M{1}.S.type,'gest')';
ix_extr = ismember(M{1}.S.type,'extr')';

Np = 5;

omega_noise = linspace(0,max_omega_range,Np);
alpha_noise = linspace(0,max_alpha_sigma,Np);

noise_pars_ixs{1} = [combvec(1:Np,1:Np)' ones(Np*Np,2)];
noise_pars_ixs{2} = [ones(Np*Np,2) combvec(1:Np,1:Np)'];

noise_pars{1} = [combvec(omega_noise,omega_noise)' zeros(Np*Np,2)];
noise_pars{2} = [zeros(Np*Np,2) combvec(alpha_noise,alpha_noise)'];


c=1;
for i=1:length(M)
    
    NP = noise_pars{noiseparset(i)};
    NP_ixs = noise_pars_ixs{noiseparset(i)};
    
    for j=1:size(NP,1)
        
        status_str = status(sprintf('model: %i/%i, params: %i/%i',i,length(M),j,size(NP,1))); %#ok<NASGU>
        
        omega_local = NP(j,1);
        omega_glob = NP(j,2);
        alpha_local = NP(j,3);
        alpha_glob = NP(j,4);      
        
        for k=1:N_batches
            
            mm = {};
            for b=1:batch_size
                m = M{i};
                
                %global and local oscillator frequency noise
                m.omega(ix_osc) = m.omega(ix_osc) ...
                    + omega_glob*(rand-1/2) ...
                    + omega_local*(rand(1,sum(ix_osc))-1/2);
                
                %global and local growth rate noise (autonomous and
                %non-autonmous):
                rn_alpha = randn;
                
                m.alpha = m.alpha ...
                    + alpha_glob*rn_alpha*(m.alpha==1) ...
                    + alpha_local*randn(size(m.alpha)).*(m.alpha==1);
                
                m.autog = m.autog ...
                    + alpha_glob*rn_alpha*ix_extr ...
                    + alpha_local*randn(size(m.autog)).*ix_extr;
                
                mm{b} = m;
            end
    
            mm = fbmod_run_models(mm,'parallel');
            
%             passed = fbmod_sanitycheck(mm,...
%                 {'precendence' 'gest_1 t_on' 'g2 t_on'},...
%                 {'precendence' 'gest_1 t_on' 'g2 t_on'});
            
            for b=1:length(mm)
                R(c,1).model = mm{b}.model_descr; %#ok<*AGROW>
                R(c).noise_pars = NP(j,:);
                R(c).omega_local = omega_local;
                R(c).omega_glob = omega_glob;
                R(c).alpha_local = alpha_local;
                R(c).alpha_glob = alpha_glob;  
                R(c).parset = noiseparset(i);
                R(c).noise_par_ixs = NP_ixs(j,:);
                [~,~,R(c).t_on,R(c).t_off] = fbmod_gestural_activation(mm{b}.X(:,ix_gest),mm{b}.t);
                c=c+1;
            end
        end
    end
end

%%
R = struct2table(R);
save([h.sim_dir mfilename '.mat'],'R');

end
