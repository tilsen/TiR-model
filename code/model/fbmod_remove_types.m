function [M] = fbmod_remove_types(M,rmtypes)

if nargin==0, return; end

rmid = ismember(M.S.type,rmtypes);

rm_vecs = {'X' 'modval' 'act_tol' 'omega' 'k' 'delta' 'floors' 'ceils' 'X0'};
rm_mats = {'phi' 'chi' 'tau' 'alpha' 'sigma' 'beta' 'map_oscr1_osc' 'map_oscr2_oscr1' 'map_exg_extr'};


for i=1:length(rm_vecs)
    M.(rm_vecs{i}) = M.(rm_vecs{i})(:,~rmid);
end

for i=1:length(rm_vecs)
    M.(rm_mats{i}) = M.(rm_mats{i})(~rmid,~rmid);
end

M.S = M.S(~rmid,:);
M.Ns = height(M.S);

end




