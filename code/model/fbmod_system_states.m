function [] = fbmod_system_states(M,X,types)

X = X(:);

M.S.x = X;
ixs = ismember(M.S.type,types);

disp(M.S(ixs,:));

end