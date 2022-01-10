function [M,U,P] = fbmod_gen_utterance_model(modelvar,sysnames,varargin)

%initialize model
allsys = [sysnames{:}];
allsys = [allsys{:}];

if numel(unique(allsys))<numel(allsys)
    fprintf('Error: non-unique system name\n');
    return;
end

pwrd_subsys = cellfun(@(c)length(c),sysnames);
Ncq = length(sysnames)+sum(pwrd_subsys);

c=length(sysnames)+1;
for i=1:length(sysnames)
    pwrd_cqsystems{i} = (c:(c+pwrd_subsys(i)-1))'; %#ok<AGROW>
    c = pwrd_cqsystems{i}(end)+1;
end

M = fbmod_init_model(allsys,varargin{:},'default_extr',false,'Ncq',Ncq);

M.S.name_model = M.S.name;
M.S.name_model(ismember(M.S.type,'gest')) = allsys';

%close all extr gates
M.X0(1,ismember(M.S.type,'exg')) = 0;

U=[];
for i=1:length(sysnames)
    [M,S] = gen_pword(modelvar,M,sysnames{i},pwrd_cqsystems{i});
    S.pwrd = i*ones(height(S),1);
    U = [U; S];     
end

U.set = U.name;
U.name = arrayfun(@(c){sprintf('sylb_%i',c)},(1:height(U))');

%% setup p-wrd gating

%associate first N cq systems with N pwrds, find last gesture
pwrds = unique(U.pwrd);
for i=1:length(pwrds)
    gest_typeids = U.gest_typeid(U.pwrd==pwrds(i));
    gest_typeids = [gest_typeids{:}];
    
    cq_gest_typeids = U.cq_gest_id(U.pwrd==pwrds(i));
    
    P(i).pwrd = i;
    P(i).cq_typeid = i;
    P(i).last_gest_typeid = cq_gest_typeids(end);    
end

P = struct2table(P);
P.cq_name = arrayfun(@(c){sprintf('cq_%i',c)},P.cq_typeid);
P.cqg_name = arrayfun(@(c){sprintf('cqg_%i',c)},P.cq_typeid);
P.last_gest_name = arrayfun(@(c){sprintf('gest_%i',c)},P.last_gest_typeid);
P.cq_id = M.getid(P.cq_name);
P.cqg_id = M.getid(P.cqg_name);

%associate cq systems with associated set of gates:
M.map_cq_cqg(P.cq_typeid,P.cq_typeid) = true(numel(P.cq_id));

%initially degate pwd cqs
M.X0(:,P.cqg_id) = 1;

%initial activation gradient:
M.X0(:,P.cq_id) = linspace(0.5,0.1,numel(P.cq_id));

sn = @(str,n)sprintf('%s_%i',str,n);
for i=1:height(U)
    
    Pix = U.pwrd(i);
    
    %connect pwd cq to within-pwd cqg, use nan threshold:
    %M = fbmod_add_action(M,P.cq_name(Pix),U.cqg_names(i),1,nan);
    cq_typeid = M.S.type_id(ismember(M.S.name,U.cqg_names{i}));
    M.map_cq1_cqg0(P.cq_typeid(Pix),cq_typeid) = 1; 
    
    %use last non-release gesture of set to supress pwd cq:
    M = fbmod_add_action(M,{sn('etmr',P.last_gest_typeid(Pix))},P.cq_name(Pix),U.cq_gest_tau(i),-1);
    
end

%make pwd cqs level=1
M.cq_level(P.cq_typeid) = 1;


end


%%
function [M,S] = gen_pword(modelvar,M,sysnames,cqsystems,varargin)

sn = @(str,n)sprintf('%s_%i',str,n);

S.name = arrayfun(@(c){sn('sylb',c)},(1:length(sysnames))');
S = struct2table(S);
for i=1:length(sysnames)
    S.gests{i} = sysnames{i};
end
S.Ng = cellfun(@(c)length(c),S.gests);
S.gest_str = cellfun(@(c)strjoin(c,' '),S.gests,'un',0);

%classify and specify initiation
%sys_ix=1;
for i=1:height(S)
    
    %S.gest_cat{i} = cellfun(@(c)regexprep(c,'_\d+',''),S.gests{i},'un',0);
    
    S.gest_cat{i} = cellfun(@(c)regexprep(c,'_{?\d+}?',''),S.gests{i},'un',0);
    
    S.cat_str(i) = cellfun(@(c){strjoin(c,' ')},S.gest_cat(i));
    S.gest_num{i} = cellfun(@(c)str2double(c),regexp(S.gest_str{i},'\d+','match'));  

    
    S.gest_id{i} = cellfun(@(c)M.S.id(ismember(M.S.name_model,c)),S.gests{i});
    S.gest_typeid{i} = arrayfun(@(c)M.S.type_id(M.S.id==c),S.gest_id{i});
    
    %sys_ix=sys_ix+S.Ng(i);
    
    Ng = S.Ng(i);
    
    %categories: C V R v c r
    S.catm{i} = [...
        ismember(S.gest_cat{i},'C');...
        ismember(S.gest_cat{i},'V');...
        ismember(S.gest_cat{i},'R');...
        ismember(S.gest_cat{i},'v');...
        ismember(S.gest_cat{i},'c');...
        ismember(S.gest_cat{i},'r')]; 
    
    S.init_osc{i} = false(Ng);
    S.init_etmr{i} = false(Ng); %src, trg
    S.init_stmr{i} = false(Ng); %src, trg
    S.term_etmr{i} = false(Ng); %src, trg
    S.term_stmr{i} = false(Ng); %src, trg
        
    %specifications of oscillator control for various syllable shapes
    if regexp(S.cat_str{i},'^C V R')
        S.init_osc{i}(1:3,1:3) = logical(eye(3));
        S.PHI{i}                = 100*[0  1 -1; 1  0  1; -1  1  0];
        S.X0_osc{i}             = [pi/16 0 -pi/16]+pi;
        
    elseif regexp(S.cat_str{i},'^C V')
        S.init_osc{i}(1:2,1:2) = logical(eye(2));
        S.PHI{i}                = 100*[0  1; 1 0];
        S.X0_osc{i}             = [pi/16 -pi/16]+pi;
        
    elseif regexp(S.cat_str{i},'^V')
        S.init_osc{i}(1,1)      = true;
        S.PHI{i}                = 0;
        S.X0_osc{i}             = pi;        
        
    elseif regexp(S.cat_str{i},'^C C V R R')
        S.init_osc{i}(1:5,1:5)  = logical(eye(5));
        !! ToDo: implement
        S.PHI{i}                = 100*[0  1 -1; 1  0  1; -1  1  0]; 
        S.X0_osc{i}             = [pi/16 -pi/16 0 pi/32 -pi/8]+pi;
    end
    
    %post-vocalic initiation options
    S.init_etmr{i} = ismember(S.gest_cat{i},{'V' 'v' 'c'})' & ismember(S.gest_cat{i},{'c' 'r'});
    
    %disallow self-initiation or post-dictive initiation:
    %S.init_etmr{i}(logical(eye(Ng))) = 0;
    S.init_etmr{i}(logical(tril(ones(Ng)))) = 0;
    
    %use latest preceding etmr to trigger (last row):
    for c=1:Ng 
        for r=Ng:-1:1
            if S.init_etmr{i}(r,c)
                S.init_etmr{i}(1:(r-1),c) = 0;
                break;
            end
        end
    end
    
    %allow same options for sensory control
    S.init_stmr{i} = S.init_etmr{i};
    
    %terminate all with etmr, stmr
    S.term_etmr{i} = logical(eye(Ng));
    S.term_stmr{i} = logical(eye(Ng));
    
   
end

%%
%by-class durations: C V R v c r
cat_durs = [0.100 0.250 0.100 0.200 0.100 0.100];

S.gest_tau = cellfun(@(c){(cat_durs*c)},S.catm);

S.term_etmr_tau = cellfun(@(c,d){c .* d},S.gest_tau,S.term_etmr);
S.term_stmr_tau = cellfun(@(c,d){c .* d},S.gest_tau,S.term_stmr);

S.init_etmr_tau = cellfun(@(c,d){c' .* d},S.gest_tau,S.init_etmr);
S.init_stmr_tau = cellfun(@(c,d){c' .* d},S.gest_tau,S.init_stmr);

%last non-release gesture (for cq supression)
S.cq_gest_id = cell2mat(cellfun(@(c,d){d(find(ismember(c,{'V' 'v' 'c'}),1,'last'))},S.gest_cat,S.gest_typeid));
S.cq_gest_tau = cellfun(@(c,d,e)c(ismember(d,e)),S.gest_tau,S.gest_typeid,num2cell(S.cq_gest_id));

S.cq_gest_tau = S.cq_gest_tau/2;

%%

Ns = height(S);

switch(modelvar)
    case 'standard'
        
        
        %first cq system of each set
        S.cq_names = arrayfun(@(c){['cq_' num2str(c(1))]},cqsystems);
        S.cqg_names = arrayfun(@(c){['cqg_' num2str(c(1))]},cqsystems);
        S.exg_names = cellfun(@(c){['exg_' num2str(c(1))]},S.gest_typeid);
        S.extr_names = cellfun(@(c){['extr_' num2str(c(1))]},S.gest_typeid);        
        
        S.cq_ids = cellfun(@(c)M.getid(c),S.cq_names);
        S.cqg_ids = cellfun(@(c)M.getid(c),S.cqg_names);
        
        %associate cq systems with associated set of gates:
        M.map_cq_cqg(cqsystems,cqsystems) = true(numel(S.cq_ids));
                               
        %initial activation gradient of cq systems
        M.X0(ismember(M.S.name,S.cq_names)) = linspace(0.5,0.01,Ns);
        
        %connect cq system to exgates:
        M = fbmod_add_action(M,S.cq_names,S.exg_names,1,1);
        
        for i=1:Ns
            
            Ng = S.Ng(i);
            gixs = S.gest_typeid{i};
            
            %oscillator initiation:         
            for a=1:Ng
                for b=1:Ng
                    if S.init_osc{i}(a,b)
                        
                        %use set-selection extrinsic system to de-gate all oscillators:
                        M = fbmod_add_action(M,S.extr_names(i),{sn('oscg',gixs(b))},0.005,1); 
                        
                        %connect oscillators to gestural gates:
                        M = fbmod_add_action(M,{sn('osc',gixs(a))},{sn('gate',gixs(b))},0,1); 
                        
                        %use internal feedback to close oscillator gates
                        M = fbmod_add_action(M,{sn('etmr',gixs(b))},{sn('oscg',gixs(b))},0.075,-1); 
                    end
                end
            end
            
            %initial oscillator phases and phase coupling:
            osc_num = S.gest_typeid{i}(any(S.init_osc{i}));
            osc_ixs = arrayfun(@(c)M.getid({sn('osc',c)}),osc_num);
            M.X0(1,osc_ixs) = S.X0_osc{i};
            M.phi(osc_ixs,osc_ixs) = S.PHI{i};            
            
            
            %etmr and stmr initiation/termination:
            for a=1:Ng
                for b=1:Ng
                    
                    %internal feedback to open gestural gate
                    if S.init_etmr{i}(a,b)                        
                        M = fbmod_add_action(M,{sn('etmr',gixs(a))},{sn('gate',gixs(b))},S.init_etmr_tau{i}(a,b),1); 
                    end
                    
                    %external feedback to open gestural gate
                    if S.init_stmr{i}(a,b)                        
                        M = fbmod_add_action(M,{sn('stmr',gixs(a))},{sn('gate',gixs(b))},S.init_stmr_tau{i}(a,b),1); 
                    end                    
                    
                    %internal feedback to close gestural gate
                    if S.term_etmr{i}(a,b)                        
                        M = fbmod_add_action(M,{sn('etmr',gixs(a))},{sn('gate',gixs(b))},S.term_etmr_tau{i}(a,b),-1); 
                    end
                    
                    %external feedback to close gestural gate
                    if S.term_stmr{i}(a,b)                        
                        M = fbmod_add_action(M,{sn('stmr',gixs(a))},{sn('gate',gixs(b))},S.term_stmr_tau{i}(a,b),-1); 
                    end                    
                    
                end
            end            
            
            %use last non-release gesture of set to supress cq:
            M = fbmod_add_action(M,{sn('etmr',S.cq_gest_id(i))},S.cq_names(i),S.cq_gest_tau(i),-1);
          
        end 
end
  

end