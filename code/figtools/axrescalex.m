function [] = axrescalex(rsvals,ax)

if nargin==0
    rsvals = [-0.01 0.01];
    ax = gca;
end

if nargin==1
    ax = gca;
end

if numel(ax)>1
    for j=1:length(ax)
        axrescalex(rsvals,ax(j));
    end
    return;
end

if numel(rsvals)==1, rsvals=rsvals*[-1 1]; end


xlims = get(ax,'xlim');
set(ax,'xlim',xlims+diff(xlims)*rsvals);

end

