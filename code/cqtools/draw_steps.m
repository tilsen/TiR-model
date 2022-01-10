function [H] = draw_steps(L,varargin)

p = inputParser;
 
default_steptype = 'sigmoid'; %sigmoid, step
default_dx = 0.001; 
default_steepness = 200; 
default_color = [0 0 0];
default_steplinewidth = 1;

addRequired(p,'L');
addOptional(p,'steptype',default_steptype);
addOptional(p,'dx',default_dx);
addOptional(p,'color',default_color);
addOptional(p,'steplinewidth',default_steplinewidth);
addOptional(p,'steepness',default_steepness);

parse(p,L,varargin{:});

extend_xx = mode(diff(L.x));

xx = min(L.x):p.Results.dx:(max(L.x)+extend_xx);

sigfcn = @(center,range)range./(1+exp(-p.Results.steepness*(xx-center)));

switch(p.Results.steptype)
    case 'sigmoid'
        Lx = [L.x; max(xx)];
        Ly = L.y;
        dy = [0; diff(Ly)];
        
        for i=2:length(Lx)-1
            Y(i,:) = sigfcn(Lx(i),dy(i));
        end
        
        yy = sum(Y);
        
end

H.lh = plot(xx,yy,'color',p.Results.color,'linew',p.Results.steplinewidth); hold on;
axis equal;

H.levelcenters = Lx(1:end-1) + diff(Lx)/2;

%draw null objects to facilitate connections
%H.levelcenterh = 

end

