function [H] = draw_cq_potential(ax,E,L,varargin)

if ~(numel(ax)==length(E)) || ~(length(E)==length(L))
    fprintf('input arguments must match in length\n'); return;
end

for i=1:height(E)
    h = draw_potential(ax(i),E{i},L{i},varargin{:});
    H(i,1) = h;
end


end

%%
function [H] = draw_potential(ax,E,L,varargin)

set(gcf,'currentaxes',ax);

H = draw_steps(L,varargin{:});
H = draw_systems_in_steps(H,E,L,varargin{:});

axis tight;
set(gca,'Visible','off');

end