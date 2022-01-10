function [] = fbmod_plot_system_matrix(M,param,varargin)

p = inputParser;

def_parent = NaN;
def_show_values = true;
def_fontsize = 18;
def_fontsize_values = 12;
def_plot_transform = {};
def_ignore_values = [0];

addRequired(p,'M');
addRequired(p,'param');
addOptional(p,'parent',def_parent);
addOptional(p,'fontsize',def_fontsize);
addOptional(p,'fontsize_values',def_fontsize_values);
addOptional(p,'show_values',def_show_values);
addOptional(p,'plot_transform',def_plot_transform);
addOptional(p,'ignore_values',def_ignore_values);

parse(p,M,param,varargin{:});

if isnan(p.Results.parent)
    figure;
else
    set(gcf,'currentaxes',p.Results.parent);
end


S = M.S;

if ischar(param)
    paramstr = param;
    param = M.(param);
else
    paramstr = '';
    param = param;
end
param = full(param);

switch(paramstr)
    case {'map_cq_cqg','map_cq1_cqg0'}
        names{1} = S.name(ismember(S.type,'cq'));
        names{2} = S.name(ismember(S.type,'cqg'));
    otherwise
        names{1} = S.name;
        names{2} = S.name;
end

if size(param,1)>1 %matrix parameter
    ix_row = ~all(ismember(param,p.Results.ignore_values) | isnan(param),2);
    ix_col = ~all(ismember(param,p.Results.ignore_values) | isnan(param),1);
    param = param(ix_row,ix_col);
    names_row = names{1}(ix_row);
    names_col = names{2}(ix_col);
else
    ix_col = ~(ismember(param,p.Results.ignore_values) | isnan(param));
    param = param(ix_col);
    names_col = names{2}(ix_col);
    names_row = {};
end

plot_param = param;
if ~isempty(p.Results.plot_transform)
    plot_param = p.Results.plot_transform(plot_param);
end

imagesc(plot_param); hold on;
matrix_gridlines([.5 .5 .5]);

clims = double(minmax(plot_param(:)'));

clims = minmax([0 clims]); %include zero;

%scale-dependent heatmap
if min(clims)==0
    cmap = stf_colormap(150,[1 1 1],[0 0 1]);
    cmap = cmap(1:end-50,:);
elseif max(clims)==0
    cmap = stf_colormap(150,[1 0 0],[1 1 1]);
    cmap = cmap(51:end,:);   
else
    cmap = stf_colormap(200,[1 0 0],[1 1 1],[0 0 1]);
    cmap = cmap(51:150,:);
end

set(gca,'XTick',1:size(param,2),'YTick',1:size(param,1),'fontsize',p.Results.fontsize,'XAxisLocation','top');
set(gca,'XTickLabel',names_col,'YTickLabel',names_row,'Ticklabelinterpreter','none','CLim',clims);
colormap(cmap);
colorbar;

if p.Results.show_values
    
    isvalid = @(x)~isnan(x) & ~ismember(x,p.Results.ignore_values);
    
    color = [0 0 0];
    
    for a=1:size(param,1)
        for b=1:size(param,2)
            if isvalid(param(a,b))
                text(b,a,num2str(param(a,b)),...
                    'fontsize',p.Results.fontsize_values,'fontweight','bold',...
                    'hori','center','verti','mid','color',color);
            end
        end
    end
end

end
