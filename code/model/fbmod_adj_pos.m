function [M] = fbmod_adj_pos(M,name,adj)

if ~iscell(name)
    names = {name};
else
    names = name;
end

for i=1:length(names)
    name = names{i};
    
    stype = M.S.type{ismember(M.S.name,name)};
    
    if isempty(stype)
        fprintf('system %s not found\n',name); return; end
    
    switch(stype)
        case {'itmr' 'etmr' 'stmr' 'osc' 'extr'}
            ixs = find(ismember(M.T.name,name));
            M.T.pos(ixs,:) =  M.T.pos(ixs,:) + adj;
            
        case 'sens'
            ixs = find(ismember(M.Y.name,name));
            M.Y.pos(ixs,:) =  M.Y.pos(ixs,:) + adj;
            
        case 'gest'
            ixs = find(ismember(M.G.name,name));
            M.G.pos(ixs,:) =  M.G.pos(ixs,:) + adj;
            
    end
end


end