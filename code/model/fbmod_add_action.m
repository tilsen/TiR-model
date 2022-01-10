function [M] = fbmod_add_action(M,varargin)

if nargin==0, return; end

%---ex:
%M=fbmod_add_action(M,src_names,targ_names,dur_thresh,valence)

%NOTE: coupling matrices: rows act on columns

%collect parameters and expand to vectors if necessary
for i=1:length(varargin)
    if ischar(varargin{i}), varargin{i} = {varargin{i}}; end
    P{i} = varargin{i}(:)';
end

%expand if any one is vector
lens = cellfun(@(c)length(c),P);
if numel(unique(lens))>2
    fprintf('error: invalid parameter inputs\n'); return;
elseif numel(unique(lens))==2
    for i=1:length(P)
        if numel(P{i})==1
            if iscell(P{i})
                P{i} = repmat(P{i},1,max(lens));
            else
                P{i} = P{i}*ones(1,max(lens));
            end
        end
    end
end

warning('off','MATLAB:table:RowsAddedExistingVars');

T = M.T;

S = M.S;
getid = M.getid;

src_ids = cellfun(@(c)getid(c),P{1});

try
    trg_ids = cellfun(@(c)getid(c),P{2});
catch
    fprintf('%s -> %s: missing system\n',P{1}{:},P{2}{:}); return;
end

dur_thresh = P{3};
valence = P{4};

for i=1:length(src_ids)
    
    Tnew = fbmod_timer_constructor(...
        S.type{src_ids(i)},...
        src_ids(i),...
        'trg',trg_ids(i),...
        'thresh',dur_thresh(i),...
        'valence',valence(i)); 
    
    T = [T; Tnew]; %#ok<*AGROW>
    
    M.tau(src_ids(i),trg_ids(i)) = dur_thresh(i);
    M.chi(src_ids(i),trg_ids(i)) = valence(i);
    
end

M.T = T;

end




