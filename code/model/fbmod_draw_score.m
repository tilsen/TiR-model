function [h] = fbmod_draw_score(G)

if ~ismember(G.Properties.VariableNames,'row')
    G.row = (1:height(G))';
end

yh = 0.90;
trng = minmax([G.t_on{:} G.t_off{:}]);
xo = 0.01*diff(trng);

h = [];

for i=1:height(G) 
    t_on = G.t_on{i};
    t_off = G.t_off{i};
    row = -G.row(i);
    for j=1:length(t_on)
        h(i).fh(j) = fill([t_on(j)*[1 1] t_off(j)*[1 1]],row+[0 yh yh 0],G.color(i,:),'facealpha',G.facealpha(i),'edgecolor','none'); hold on;
        h(i).th(j) = text(t_on(j)+xo,row+yh/2,G.label{i},'fontsize',G.fontsize(i),'interp',G.interpreter{i});
    end  
end

set(gca,'YColor','none','YTick',[],'Box','off');
axis tight; axrescaley(0.05);

end