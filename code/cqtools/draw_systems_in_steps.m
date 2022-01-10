function [H] = draw_systems_in_steps(H,E,L,varargin)

p = inputParser;
p.KeepUnmatched = true;

def_include_operations = true;
def_operation_arrowprops = {'length',4,'tipangle',30};
 
addRequired(p,'H');
addRequired(p,'E');
addRequired(p,'L');
addOptional(p,'include_operations',def_include_operations);
addOptional(p,'operation_arrowprops',def_operation_arrowprops);

parse(p,H,E,L,varargin{:});

th = linspace(0,2*pi,100);

for i=1:height(L)
    ixs = L.ix_sys(i);
    
    if iscell(ixs), ixs=ixs{1}; end
    
    if isnan(ixs), continue; end
    
    xc = H.levelcenters(i);
    
    for j=1:length(ixs)
        
        ee = E(ixs(j),:);
        
        xx = xc + ee.radius*cos(th);
        yy = L.y(i) + ee.radius + ee.radius*sin(th); 
        
        H.sys_fh(i,j) = fill(xx,yy,ee.color,'facealpha',ee.facealpha,'edgecolor','k');
        
        switch(ee.textloc{:})
            case 'inside'
                
                H.sys_th(i,j) = text(xc,L.y(i)+ee.radius,ee.label,'interpreter',ee.interpreter{:},...
                    'hori','center','verti','mid','fontsize',ee.fontsize);
                
            case 'above'
                 
                H.sys_th(i,j) = text(xc,L.y(i)+2*ee.radius,ee.label,'interpreter',ee.interpreter{:},...
                    'hori','center','verti','bot','fontsize',ee.fontsize);               
                
        end
        
        if p.Results.include_operations
            switch(ee.op)
                case 1
                    PP = [xc max(yy); nan nan; L.x(i+1) L.y(i+1)+2*ee.radius];
                    PP(2,:) = [PP(1,1) PP(3,2)];
                    H.op_arh{i,j} = draw_arc(PP,'arrowprops',p.Results.operation_arrowprops);
                    
                case 0
                    H.op_arh{i,j} = [];
            
                case -1
                    dy = 0.05*diff(L.y(1:2));
                    PP = [H.levelcenters(end) L.y(end); nan(2); xc L.y(i)-dy];
                    PP(2,:) = [H.levelcenters(end-1) L.y(i)]; 
                    PP(3,:) = [PP(end,1) L.y(i)] - [-diff(L.x(1:2)) diff(L.y(1:2))]/2;
                    H.op_arh{i,j} = draw_arc(PP,'arrowprops',p.Results.operation_arrowprops);                    
            end
        else
            H.op_arh{i,j} = [];
        end
    end
    
end

end

