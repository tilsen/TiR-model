function [L] = gen_cq_levels(E,varargin)

%loop over table rows
if height(E)>1
    for i=1:height(E)
        L(i,1) = gen_cq_levels(E(i,:),varargin{:});
    end
    return;
end
    
p = inputParser;

def_levelspacingmode = 'equal';
def_actspacingmode = 'mirror';
def_includeselectionlevel = true;
def_includegroundlevel = true;
def_groundvalue = 0;
def_selectionvalue = 1;

addRequired(p,'E',@(x)istable(E));
addOptional(p,'levelspacingmode',def_levelspacingmode);
addOptional(p,'actspacingmode',def_actspacingmode);
addOptional(p,'includeselectionlevel',def_includeselectionlevel);
addOptional(p,'includegroundlevel',def_includegroundlevel);
addOptional(p,'groundvalue',def_groundvalue);
addOptional(p,'selectionvalue',def_selectionvalue);

parse(p,E,varargin{:});

level_vals = E.X_ordered;

if p.Results.includegroundlevel
    level_vals = unique([p.Results.groundvalue level_vals]);
end

if p.Results.includeselectionlevel
    level_vals = unique([p.Results.selectionvalue level_vals]);
end

N_levels = length(level_vals);

switch(p.Results.levelspacingmode)
    case 'equal'
        L.y = linspace(min(level_vals),max(level_vals),N_levels)';
        
    case 'data'
        L.y = level_vals';
end

switch(p.Results.actspacingmode)
    case 'mirror'
        L.x = L.y;
        
    case 'data'
        L.x = level_vals';
end


L = struct2table(L);

highest_is_sel = E.X_ordered(end)>=p.Results.selectionvalue;
lowest_is_ground = E.X_ordered(1)<=p.Results.groundvalue;

if height(L)==2 && lowest_is_ground %all systems grounded
    L.ix_sys = {E.ix_ordered; nan};
   
elseif height(L)==2 && highest_is_sel %all systems grounded
    L.ix_sys = {nan; E.ix_ordered};
    
elseif highest_is_sel && lowest_is_ground
    L.ix_sys = E.ix_ordered';
    
elseif highest_is_sel && ~lowest_is_ground
    L.ix_sys = [nan; E.ix_ordered'];
    
elseif ~highest_is_sel && lowest_is_ground
    L.ix_sys = [E.ix_ordered'; nan];
    
elseif ~highest_is_sel && ~lowest_is_ground
    L.ix_sys = [nan; E.ix_ordered'; nan];
    
end

L = {L};

end