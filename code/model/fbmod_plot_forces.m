function [ph,Fh] = fbmod_plot_forces(T,M,h)

if isempty(T)
    T = M.T(~cellfun('isempty',M.T.F),:); 
end

for i=1:height(T)
    if isempty(T.F{i}), continue; end
    for j=1:size(T.F{i},2)
        ph(i,j) = plot(M.t,T.F{i}(:,j),'linew',2,'color',T.color(i,:)); hold on;
    end
end

Fh = fbmod_label_actions(T,M);
ylim([-1.5 1.5]);
set(gca,'ytick',[-1 1],'YTickLabel',h.fontsizer({'-','+'}',h.fs(2)),'XTickLabel',[],...
    'Box','off','ygrid','on');


end