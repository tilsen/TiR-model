function [M] = fbmod_01(M)

catch_events = false;
catch_times = [];

M.limit = @(x,floors,ceils)min([max([x; floors]); ceils]);
M.actfcn = @(x,y)sqrt((x-y).^2)<(M.act_tol');

M.fXX = @(x)(M.chi.*M.actfcn(x',M.tau));
M.FXX = @(x)nansum(M.fXX(x),1); %#ok<*NANSUM>

M.phiXX = @(x)sin(x'-x).*M.phi;
M.PHIXX = @(x)nansum(M.phiXX(x),1);

%put some variables in workspace
%ff = fieldnames(M);
ff = {'t','dt','X','Nt','Ns','delta'};
for i=1:length(ff)
    eval([ff{i} '=M.(ff{i});']);
end

M.F_XX = zeros(Nt,Ns,Ns);
M.PHI_XX = zeros(Nt,Ns,Ns);

%initial conditions
X(1,:) = M.X0;
notisinfmodval = ~isinf(M.modval);
X(1,notisinfmodval) = mod(X(1,notisinfmodval),M.modval(notisinfmodval));

%auxiliary indices
isosc   = ismember(M.S.type,'osc')'; 
issens  = ismember(M.S.type,'sens')';
isgest  = ismember(M.S.type,'gest')';
% isextr  = ismember(M.S.type,'extr')';
% isstmr  = ismember(M.S.type,'stmr')';
% isetmr  = ismember(M.S.type,'etmr')';
% isitmr  = ismember(M.S.type,'itmr')';
iscq    = ismember(M.S.type,'cq')'; 
iscqg   = ismember(M.S.type,'cqg')'; ixcqg = find(iscqg);

M.cq_sel = false(1,sum(iscq));

tic

for i=2:Nt
    
    if any(abs(t(i)-catch_times)<dt)
        %for debugging at specific time, set catch_times and put breakpoint here
        %fbmod_system_states(M,X(i-1,:),{'cq' 'cqg'});
        1;
    end
   
    xx = X(i-1,:);
    
    %competitive queuing: 
    M.cq_sel(xx(iscq)>=1) = true; %flag to indicate if cq has been selected
    
    %cond1: any cq in set is currently selected
    cond1 = any((xx(iscq)'.*M.map_cq_cqg)>=1); 
    
    %cond2: all cqs in set previously selected
    num_cq_prev_sel = sum(M.cq_sel'.*M.map_cq_cqg);
    cond2 = (num_cq_prev_sel~=0) & (num_cq_prev_sel==sum(M.map_cq_cqg)); 
    
    %cond3: higher-level cq system exists but is not selected:
    cq_cqg = M.map_cq1_cqg0; 
    higher_lev_cq_sel = (xx(iscq).*(M.cq_level==1)>=1);
    cond3 = (sum(higher_lev_cq_sel'.*cq_cqg)==0) & ~(M.cq_level==1);
    
    %open gates by default:
    xx(iscqg) = 1;
    
    %close gates if any of three conditions are met:
    xx(ixcqg(cond1 | cond2 | cond3)) = 0;
         
    %oscillator action forces depend on phase only when amplitude > 
    xxi = xx;
    osc_oscr1_thresh = (xx>=M.r1_thresh)*M.map_oscr1_osc;
    xxi(isosc & ~osc_oscr1_thresh) = nan;
    
    %force records:
    M.F_XX(i,:,:) = M.fXX(xxi); %#ok<*UNRCH>
    M.PHI_XX(i,:,:) = M.phiXX(xx);
    
    %catch force events on activation control for debugging
    if catch_events
        
        forces = squeeze(M.F_XX(i,:,:))~=0; 
        if any(forces,'all')
            [ixr,ixc] = find(forces);
            for k=1:length(ixr)
                fprintf('t=%1.3f\t%s -> %s\t%1.1f\n',...
                    t(i),M.S.name{ixr(k)},M.S.name{ixc(k)},M.F_XX(i,ixr(k),ixc(k)));
            end 
            1;
        end
    end
    
    %oscillator amplitude dynamics:
    dosc_amp = (xx*M.map_oscr2_oscr1) + (-2*sqrt(M.k).*xx - M.k.*(xx*(M.map_oscr2_oscr1')-xx*M.beta));
       
    %Eq. (1):
    dX =  M.omega.*(xx*M.map_exg_extr) + ...                                           %gated autonomous growth external systems
                    dosc_amp + ...                                                     %oscillator amplitude dynamics
          M.omega.*(xx*M.map_oscr1_osc).^2 + M.PHIXX(xx) + ...                         %oscillator phase dynamics
          xx*M.alpha + (1-xx)*M.sigma + ...                                            %gestural gating forces, feedback integration, and supression
          (1/dt).*M.FXX(xxi);                                                          %action-forces      
    
    X(i,:) = xx + dX*dt;
    
    %sensory systems
    delayed_tix = max(i-delta(issens),1);
    X(i,issens) = diag(X(delayed_tix,isgest))';
    
    X(i,notisinfmodval) = mod(X(i,notisinfmodval),M.modval(notisinfmodval));
    X(i,:) = M.limit(X(i,:),M.floors,M.ceils);
        
end

M.elapsed_time = toc;

%%

M.X = X;

end