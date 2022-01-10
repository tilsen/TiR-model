function [M] = fbmod_remove_systems(M,rmtypes,rmsys)

if nargin==0, return; end
if nargin==2, rmsys = {}; end

rmid = ismember(M.S.type,rmtypes);
rmid = rmid | ismember(M.S.name,rmsys);

rm_vecs = {'X' 'X0' 'modval' 'act_tol' 'omega' 'k' 'delta' 'floors' 'ceils'};
rm_mats = {'phi' 'chi' 'tau' 'alpha' 'sigma' 'beta' 'map_oscr1_osc' 'map_oscr2_oscr1' 'map_exg_extr'};


for i=1:length(rm_vecs)
    M.(rm_vecs{i}) = M.(rm_vecs{i})(:,~rmid);
end

for i=1:length(rm_mats)
    M.(rm_mats{i}) = M.(rm_mats{i})(~rmid,~rmid);
end

M.S = M.S(~rmid,:);
M.Ns = height(M.S);
M.S.id = (1:M.Ns)';

%reconstruct timers
[rix,cix] = find(M.chi~=0);
T = [];
for i=1:length(rix)
    Tn = fbmod_timer_constructor(...
        M.S.type{rix(i)},...
        rix(i),...
        'trg',cix(i),...
        'thresh',M.tau(rix(i),cix(i)),...
        'valence',M.chi(rix(i),cix(i)));
    
    T = [T; Tn]; %#ok<*AGROW>
end

M.T = sortrows(T,{'src'});


end




