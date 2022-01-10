function [D,E] = gen_cq_config(t,X,names,varargin)

p = inputParser;

%force time=rows
t = t(:);

def_sel_thresh = 1;
def_output_epochs = 'anyselected'; %all

addRequired(p,'t',@(x)isvector(x));
addRequired(p,'X',@(x)ismatrix(x) & length(x)==length(t));
addRequired(p,'names',@(c)iscell(c));
addOptional(p,'sel_thresh',def_sel_thresh);
addOptional(p,'output_epochs',def_output_epochs);

parse(p,t,X,names,varargin{:});

%system is selected
is_sel = X>=p.Results.sel_thresh;

%any change in selection from previous frame
change_sel = [1; find(any(diff(is_sel),2))+1];

for i=1:length(change_sel)
    D(i).t = t(change_sel(i));
    D(i).any_sel = any(is_sel(change_sel(i),:));
    D(i).is_sel = is_sel(change_sel(i),:);  
    D(i).ix_sel = find(D(i).is_sel);
    D(i).X = X(change_sel(i),:);
    [D(i).X_ordered,D(i).ix_ordered] = sort(D(i).X);
    D(i).sel_ordered = D(i).X_ordered>=p.Results.sel_thresh;
    D(i).names = names(:)';
    D(i).names_ordered = names(D(i).ix_ordered)';
    
    %operation information
    D(i).op = zeros(1,numel(D(i).names));
    if ~any(D(i).any_sel)
        if i<length(change_sel)
            D(i).op = ones(1,numel(D(i).names));
        end
        if i>1
            D(i).op(D(i-1).is_sel) = -1;
        end
    end
    
end

D = struct2table(D);

n_rows = size(X,2);
for i=1:height(D)
    
    ff = D.Properties.VariableNames;
    ff = ff(~contains(ff,'_ordered'));
    ee = table2struct(D(i,ff));
    
    for j=1:length(ff)
        if size(ee.(ff{j}),2)>size(ee.(ff{j}),1)
            ee.(ff{j}) = ee.(ff{j})';
        else
            ee.(ff{j}) = repmat(ee.(ff{j}),n_rows,1);
        end
    end
    
    if length(ee.t)==1
        E{i,1} = struct2table(ee,'AsArray',true);
    else
        E{i,1} = struct2table(ee);
    end
    
end

switch(p.Results.output_epochs)
    case 'anyselected'
        E = E(D.any_sel);
        D = D(D.any_sel,:);
end

end