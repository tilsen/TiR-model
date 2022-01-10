function [H] = fbmod_draw_connections(H,C,varargin)

p = inputParser;

def_arrow_props = {'length',6,'tipangle',30};
def_linewidth = 1;

addRequired(p,'H');
addRequired(p,'C');
addOptional(p,'arrow_props',def_arrow_props);
addOptional(p,'linewidth',def_linewidth);

parse(p,H,C,varargin{:});

%draw
warning('off','arrow:limits');
for i=1:height(C)
    if any(isnan(C.conn{i}(:))), continue; end
    H(end+1).name = C.name{i};
    H(end).type = 'connection';
    
    H(end).handles = draw_connection(...
        C.conn{i}',...
        'linewidth',C.linewidth(i),...
        'conntype',C.conntype{i},...
        'arrowprops',p.Results.arrow_props);
    
    H(end).hvec = H(end).handles;
    H(end).conntype = C.conntype{i};
    
    H(end).obj1 = C.obj1(i);
    H(end).obj2 = C.obj2(i);
    H(end).obj1_vis = H(end).obj1.Visible;
    H(end).obj2_vis = H(end).obj2.Visible;
    
    if ismember('off',[H(end).obj1_vis  H(end).obj2_vis]) 
        set(H(end).hvec,'Visible','off'); 
    end
    
end

end


