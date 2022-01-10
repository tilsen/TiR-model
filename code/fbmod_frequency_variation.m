function [] = fbmod_frequency_variation()

dbstop if error;
h = fbmod_helpers();

%% determine frequencies

frq_range = 8*(2*pi); %range of frequencies (in rad/s)
lambda_resp = 5; %sensitivity of frequency to lambda
base_frq = 6*(2*pi); %baseline frequency

%function related frequency to lambda
frq_eff = @(lambda)base_frq - (frq_range/2)*tanh(lambda_resp*(lambda-1/2));


lambdas = linspace(0,1,21);
omegas = frq_eff(lambdas); %simulation frequencies


%% specify and run models
m = fbmod_models('CV model');
osc_ix = ismember(m.S.type,'osc');

for i=1:length(lambdas)
    m.omega(osc_ix) = omegas(i);
    M{i} = m;    
    R(i,1).model = sprintf('freqvar_%02i',i);
    R(i).lambda = lambdas(i);
    R(i).omega = omegas(i);
end
M = fbmod_run_models(M);

%% extract gestural timing
gest_ix = ismember(M{1}.S.type,'gest');
for i=1:length(M)
    [~,~,R(i).t_on,R(i).t_off] = fbmod_gestural_activation(M{i}.X(:,gest_ix),M{i}.t);
end
R = struct2table(R);
R.t_on = cell2mat(R.t_on);
R.t_off = cell2mat(R.t_off);
R.deltas = diff(R.t_on,[],2);

M = fbmod_prep_schema(M([1 floor(length(M)/2) end]));

%%

save([h.sim_dir mfilename '.mat'],'R','M');


end

