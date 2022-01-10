function [M,S] = fbmod_gen_pword_model(modelvar,sysnames,varargin)


sn = @(str,n)sprintf('%s_%i',str,n);

S.name = arrayfun(@(c){sn('sylb',c)},(1:length(sysnames))');
S = struct2table(S);
for i=1:length(sysnames)
    S.gests{i} = sysnames{i};
end
S.Ng = cellfun(@(c)length(c),S.gests);
S.gest_str = cellfun(@(c)strjoin(c,' '),S.gests,'un',0);

%classify and specify initiation
sys_ix=1;
for i=1:height(S)
    
    S.gest_cat{i} = cellfun(@(c)regexprep(c,'_\d+',''),S.gests{i},'un',0);
    S.cat_str(i) = cellfun(@(c){strjoin(c,' ')},S.gest_cat(i));
    S.gest_num{i} = cellfun(@(c)str2double(c),regexp(S.gest_str{i},'\d+','match'));    
    S.gest_ix{i} = sys_ix:(sys_ix+S.Ng(i)-1);
    
    sys_ix=sys_ix+S.Ng(i);
    
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
        
    %hand-specify oscillator control for various syllable shapes
    if regexp(S.cat_str{i},'^C V R')
        S.init_osc{i}(1:3,1:3) = logical(eye(3));
        S.PHI{i}                = 100*[0  1 -1; 1  0  1; -1  1  0];
        S.X0_osc{i}             = [pi/16 0 -pi/16]+pi;
        
    elseif regexp(S.cat_str{i},'^C C V R R')
        S.init_osc{i}(1:5,1:5)  = logical(eye(5));
        !! ToDo: implement
        S.PHI{i}                = 100*[0  1 -1; 1  0  1; -1  1  0]; 
        S.X0_osc{i}             = [pi/16 -pi/16 0 pi/32 -pi/8]+pi;
    end
    
    %post-vocalic initiation options
    S.init_etmr{i} = ismember(S.gest_cat{i},{'V' 'v' 'c'})' & ismember(S.gest_cat{i},{'c' 'r'});
    
    %disallow self-initiation:
    S.init_etmr{i}(logical(eye(Ng))) = 0;
    
    %use latest etmr to trigger (last row):
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
S.cq_gest_ix = cell2mat(cellfun(@(c,d){d(find(ismember(c,{'V' 'v' 'c'}),1,'last'))},S.gest_cat,S.gest_ix));
S.cq_gest_tau = cellfun(@(c,d,e)c(ismember(d,e)),S.gest_tau,S.gest_ix,num2cell(S.cq_gest_ix));

S.cq_gest_tau = S.cq_gest_tau/2;

%%
sys = [sysnames{:}];
M = fbmod_init_model(sys,varargin{:},'default_extr',false);

Ns = height(S);

switch(modelvar)
    case 'standard'
        
        %close all extr gates except first
        %ix_exg = setdiff(find(ismember(M.S.type,'exg')),find(ismember(M.S.type,'exg'),1,'first'));
        %M.X0(1,ix_exg) = 0;
        
        %close all extr gates
        M.X0(1,ismember(M.S.type,'exg')) = 0;
        
        %first cq system of each set
        S.cq_names = cellfun(@(c){['cq_' num2str(c(1))]},S.gest_ix);
        S.cqg_names = cellfun(@(c){['cqg_' num2str(c(1))]},S.gest_ix);
        S.exg_names = cellfun(@(c){['exg_' num2str(c(1))]},S.gest_ix);
        S.extr_names = cellfun(@(c){['extr_' num2str(c(1))]},S.gest_ix);        
        
        S.cq_ids = cellfun(@(c)M.getid(c),S.cq_names);
        S.cqg_ids = cellfun(@(c)M.getid(c),S.cqg_names);
        
        %associate cq systems with associated set of gates:
        M.map_cq_cqg(S.cq_ids,S.cqg_ids) = true(numel(S.cq_ids));
        
        %set growth rates of ununsed cq systems to 0?
        %M.omega(ismember(M.S.type,'cq') & ~ismember(M.S.name,S.cq_names)) = 0;
                        
        %initially degate cq systems (2nd extr in set)
        M.X0(ismember(M.S.name,S.cqg_names)) = 1;
        
        %initial activation gradient of cq systems
        M.X0(ismember(M.S.name,S.cq_names)) = linspace(0.5,0.01,Ns);
        
        %connect cq system to exgates:
        M = fbmod_add_action(M,S.cq_names,S.exg_names,1,1);
        
        for i=1:Ns
            
            Ng = S.Ng(i);
            gixs = S.gest_ix{i};
            
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
            osc_num = S.gest_ix{i}(any(S.init_osc{i}));
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
            M = fbmod_add_action(M,{sn('etmr',S.cq_gest_ix(i))},S.cq_names(i),S.cq_gest_tau(i),-1);
          
        end 
end
  

end