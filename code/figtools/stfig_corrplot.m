function [H] = stfig_corrplot(X,ax,varargin)

p = inputParser;

default_rows = 'pairwise';
default_ax = '';

addRequired(p,'X',@(x)ismatrix(x) & size(x,1)>=size(x,2));
addOptional(p,'ax',default_ax,@(x)all(ishandle(x),'all'));
addParameter(p,'rows',default_rows);

parse(p,X,ax,varargin{:});

Nvar = size(X,2);

switch(isempty(p.Results.ax))
    case 1
        H.ax = stf([Nvar Nvar],[0.05 0.05 0.01 0.05],[0.01 0.01],'handlearray','matrix');
    otherwise
        if size(ax,1)~=size(ax,2)
            H.ax = reshape(ax,Nvar,[]);
        else
            H.ax = p.Results.ax;
        end
end

C = corr(X,'rows',p.Results.rows);

for a=1:Nvar
    for b=1:Nvar
        set(gcf,'currentaxes',H.ax(b,a));
        if a==b
            H.histh(a) = histogram(X(:,a)); hold on;
        else
            H.ph(a,b) = scatter(X(:,a),X(:,b),'ko'); hold on;
        end
    end
end

end