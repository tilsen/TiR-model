function [th] = stfig_panlab(ax,labs,varargin)

defaultlocation = 'northwest';
defaultverticalalignment = 'bottom';
defaulthorizontalalignment = 'right';
defaultfontsize = 30;
defaultfontweight = 'bold';
defaultyoffset = 0;
defaultxoffset = 0.01;
defaultpositionmode = 'dataunits';

p = inputParser;
addRequired(p,'ax',@(x)all(ishandle(x)));
addRequired(p,'labs');
addParameter(p,'location',defaultlocation);
addParameter(p,'horizontalalignment',defaulthorizontalalignment);
addParameter(p,'verticalalignment',defaultverticalalignment);
addParameter(p,'fontsize',defaultfontsize);
addParameter(p,'fontweight',defaultfontweight);
addParameter(p,'yoffset',defaultyoffset);
addParameter(p,'xoffset',defaultxoffset);
addParameter(p,'positionmode',defaultpositionmode);

parse(p,ax,labs,varargin{:});
r = p.Results;

if ischar(labs)
    labs = {labs};
end

%expand yoffset and xoffset 
if numel(r.yoffset)~=length(ax)
    r.yoffset = r.yoffset(:)';
    r.yoffset = repmat(r.yoffset,length(ax),1);
end

if numel(r.xoffset)~=length(ax)
    r.xoffset = r.xoffset(:)';
    r.xoffset = repmat(r.xoffset,length(ax),1);
end

switch(p.Results.positionmode)
    case 'outerposition'
        axbak = stbgax;
end

for i=1:length(ax)
    if isempty(labs{i}), continue; end
   
    switch(p.Results.positionmode)
        case 'outerposition'
            xlims = ax(i).OuterPosition([1 3])*[1 0; 1 1]';
            ylims = ax(i).OuterPosition([2 4])*[1 0; 1 1]';
            textax = axbak;
            scaleax = axbak;
    
        otherwise
            xlims = ax(i).XLim;
            ylims = ax(i).YLim;
            textax = ax(i);
            scaleax = ax(i);
    end
            
    switch(r.location)
        case 'northwest'
            xpos = xlims(1);
            ypos = ylims(2);
        case 'northeast'
            xpos = xlims(2);
            ypos = ylims(2);
        case 'north'
            xpos = mean(xlims);
            ypos = ylims(2);
            r.horizontalalignment = 'center';
        case 'southwest'
            xpos = xlims(1);
            ypos = ylims(1);
        case 'southeast'
            xpos = xlims(2);
            ypos = ylims(1);
        case 'south'
            xpos = mean(xlim);
            ypos = ylims(1);
            r.horizontalalignment = 'center';
    end

    
    th(i) = text(xpos,ypos,labs{i},...
        'fontsize',r.fontsize,'hori',r.horizontalalignment,'verti',r.verticalalignment,'parent',textax,'fontweight',r.fontweight);
    
    th(i).Position(2) = th(i).Position(2) + r.yoffset(i)*diff(scaleax.YLim);
    th(i).Position(1) = th(i).Position(1) + r.xoffset(i)*diff(scaleax.XLim);
end
        

end