function [H,C,h] = fbmod_draw_schema(M,h,varargin)

p = inputParser;

def_linewidth = 2;
def_reduce_lines = [0.015 0.015];
def_format_axes = true;
def_arrow_props = {'length',6,'tipangle',30};

addRequired(p,'M');
addRequired(p,'h');
addOptional(p,'linewidth',def_linewidth);
addOptional(p,'reduce_lines',def_reduce_lines);
addOptional(p,'format_axes',def_format_axes);
addOptional(p,'arrow_props',def_arrow_props);

parse(p,M,h,varargin{:});

H = fbmod_draw_gestures([],M.G);

if p.Results.format_axes
    if strcmp(get(gca,'XLimMode'),'auto')
        axis tight; axis 'equal';
    end
    set(gca,'Visible','off');
end

H = fbmod_draw_timers(H,M.T);
H = fbmod_draw_sens(H,M.Y);

C = fbmod_connect_objs(H,M,...
    'reduce_lines',p.Results.reduce_lines,...
    'linewidth',p.Results.linewidth);

H = fbmod_draw_connections(H,C,'arrow_props',p.Results.arrow_props,'linewidth',p.Results.linewidth);

H = struct2table(H);

h.getnode = @(name)table2struct(H(ismember(H.name,name),:));
h.getconn = @(name)table2struct(H(ismember(H.name,name),:));

h.getnodeh = @(name)H.hvec{ismember(H.name,name)};
h.getconnh = @(name)H.hvec{ismember(H.name,name)};


end

