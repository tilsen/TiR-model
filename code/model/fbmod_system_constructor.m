function [M] = fbmod_system_constructor(M,systype,typeid,varargin)

if ~isfield(M,'S')
    warning('off','MATLAB:table:PreallocateCharWarning');
    M.S = table('size',[0 6],...
        'variablenames',{'name','type','type_id','id','x0','descr'},...
        'variabletypes',{'char','char','int16','int16','double','char'});
end

id = size(M.X,2)+1;

M.X(:,end+1) = zeros(size(M.X,1),1);

S.name = {sprintf('%s_%i',systype,typeid)};
S.type = {systype};
S.type_id = typeid;
S.id = id;
S.x0 = 0;
S.descr = {''};
if ~isempty(varargin)
    S.descr = {varargin{1}};
end
M.S = [M.S; struct2table(S)];

end