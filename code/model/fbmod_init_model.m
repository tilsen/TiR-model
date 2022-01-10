function [M] = fbmod_init_model(names,varargin)


M.name = names;
M.Ng = length(names); Ng = M.Ng;

for i=1:2:length(varargin)
    M.(varargin{i}) = varargin{i+1};
end

def_params = {...
    'dt' 0.001;
    'dur' 1; 
    'color' lines(Ng)'; 
    'Nextr' Ng; 
    'Nosc' Ng; 
    'osc_freq' 5; 
    'default_extr' true;
    'Ncq' 0;
    'sigma_autog' 0;
    'sigma_omega' 0;
    'sigma_stmr' 0;
    'sigma_etmr' 0;
    'sigma_itmr' 0;
    };

for i=1:size(def_params,1)
    if ~isfield(M,def_params{i,1})
        M.(def_params{i,1}) = def_params{i,2};
    end
end

Nosc = M.Nosc;
Nextr = M.Nextr;
Ncq = M.Ncq;

M.t = (0:M.dt:M.dur)';
M.Nt = length(M.t); Nt = M.Nt;

M.X = zeros(Nt,0);

for i=1:Ng,  M = fbmod_system_constructor(M,'gate',i); end
for i=1:Ng,  M = fbmod_system_constructor(M,'gest',i); end
for i=1:Ng,  M = fbmod_system_constructor(M,'itmr',i); end
for i=1:Ng,  M = fbmod_system_constructor(M,'etmr',i); end
for i=1:Ng,  M = fbmod_system_constructor(M,'stmr',i); end
for i=1:Ng,     M = fbmod_system_constructor(M,'sens',i); end
for i=1:Nextr,  M = fbmod_system_constructor(M,'extr',i); end
for i=1:Nextr,  M = fbmod_system_constructor(M,'exg',i); end
for i=1:Ncq,    M = fbmod_system_constructor(M,'cq',i); end
for i=1:Ncq,    M = fbmod_system_constructor(M,'cqg',i); end
for i=1:Nosc,  M = fbmod_system_constructor(M,'osc',i); end
for i=1:Nosc,  M = fbmod_system_constructor(M,'oscr1',i); end
for i=1:Nosc,  M = fbmod_system_constructor(M,'oscr2',i); end
for i=1:Nosc,  M = fbmod_system_constructor(M,'oscg',i); end

%
M.S.x0(ismember(M.S.type,'exg')) = 1;

%utility functions
M.Ns = height(M.S);
S = M.S;
Ns = M.Ns;
dt = M.dt;

%utility functions
M.getid = @(name)S.id(ismember(S.name,name)); %get system id
M.typeixs = @(systype)ismember(S.type,systype);
getid = M.getid;
typeixs = M.typeixs;
sysn = @(s,i)sprintf('%s_%i',s,i);

%modulus value for oscillator phase
M.modval = inf*ones(1,Ns);
M.modval(ismember(M.S.type,'osc')) = 2*pi;

%
M.act_tol = (2*dt)*ones(1,Ns);
M.act_tol(ismember(M.S.type,'osc')) = 2/(2*pi*max(M.osc_freq));

%autnomous system growth rates/oscillator phase velocity
M.omega = zeros(1,Ns);
M.omega(typeixs('extr')) = 1;
M.omega(typeixs('cq')) = 10;
M.omega(typeixs('osc')) = M.osc_freq*2*pi;

%oscillator amplitude stiffness
M.k = zeros(1,Ns);
M.k(typeixs('oscr2')) = 5000;

%relative phase coupling
M.phi = sparse(zeros(Ns));

%short-timescale force matrices
M.chi = sparse(zeros(Ns));
M.tau = sparse(nan(Ns));

%interactions
alpha = 1/dt;
M.alpha = sparse(zeros(Ns));
M.beta = sparse(zeros(Ns));
M.map_oscr1_osc = sparse(zeros(Ns));
M.map_oscr2_oscr1 = sparse(zeros(Ns));

M.map_exg_extr = sparse(zeros(Ns));

%gates --> systems (rows act on columns)
for i=1:Nextr,  M.map_exg_extr(getid(sysn('exg',i)),getid(sysn('extr',i))) = 1; end
for i=1:Ncq,    M.map_exg_extr(getid(sysn('cqg',i)),getid(sysn('cq',i))) = 1; end

%cq gate mutual exclusion relations
M.map_cq_cqg = sparse(false(Ncq));
M.cq_level = zeros(1,Ncq);

%level 1 cq to level 0 cqg
M.map_cq1_cqg0 = sparse(false(Ncq));

%amplitude threshold for phase-based triggering
M.r1_thresh = 0.50;

%gates --> gestures (rows act on columns)
for i=1:Ng, M.alpha(getid(sysn('gate',i)),getid(sysn('gest',i))) = alpha; end

%oscillator gates --> oscillator amplitudes/phases
for i=1:Ng 
    M.beta(getid(sysn('oscg',i)),getid(sysn('oscr2',i))) = 1; 
    M.beta(getid(sysn('oscg',i)),getid(sysn('osc',i))) = 1; 
end

%oscillator amplitudes --> phases
for i=1:Ng 
    M.map_oscr1_osc(getid(sysn('oscr1',i)),getid(sysn('osc',i))) = 1; 
end

%oscillator amplitude value and derivative
for i=1:Ng 
    M.map_oscr2_oscr1(getid(sysn('oscr2',i)),getid(sysn('oscr1',i))) = 1; 
end
   
%gestures/sensory systems --> timers
for i=1:Ng
    M.alpha(getid(sysn('gest',i)),getid(sysn('itmr',i))) = 1; 
    M.alpha(getid(sysn('gest',i)),getid(sysn('etmr',i))) = 1; 
    M.alpha(getid(sysn('sens',i)),getid(sysn('stmr',i))) = 1; 
end
   
%supression rates
sigma = -1/dt;
M.sigma = zeros(Ns,Ns);
for i=1:Ng, M.sigma(getid(sysn('gate',i)),getid(sysn('gest',i))) = sigma; end
for i=1:Ng, M.sigma(getid(sysn('gate',i)),getid(sysn('itmr',i))) = sigma; end
for i=1:Ng, M.sigma(getid(sysn('gate',i)),getid(sysn('etmr',i))) = sigma; end
for i=1:Ng, M.sigma(getid(sysn('gate',i)),getid(sysn('stmr',i))) = sigma; end

%motor-sensory delay
M.delta_YX = 0.040;
M.delta_YX_ix = uint16(round(M.delta_YX/M.dt));
M.delta = uint16(zeros(1,Ns));
M.delta(ismember(M.S.type,'sens')) = M.delta_YX_ix;

%table for timers
M.T = fbmod_timer_constructor(...
    'extr',M.getid('extr_1'),...
    'trg',M.getid('gate_1'),...
    'thresh',M.delta_YX,'valence',1,'auto',1);

if M.default_extr
    %external trigger: ix, time (s), type (-1 or 1)
    M.tau(getid('extr_1'),getid('gate_1')) = M.delta_YX;
    M.chi(getid('extr_1'),getid('gate_1')) = 1;
    
else
    M.T = M.T([],:);
    
end

M.floors = zeros*ones(1,Ns);
M.floors(ismember(S.type,{'osc'})) = -2*pi;
M.floors(ismember(S.type,{'oscr2'})) = -inf;

M.ceils = inf*ones(1,Ns);
M.ceils(ismember(S.type,{'cqg' 'cq' 'exg' 'oscg' 'gate' 'gest' 'sens'})) = 1;

M.X0(1,:) = M.S.x0;



end


