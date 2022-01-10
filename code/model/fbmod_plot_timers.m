function [ph] = fbmod_plot_timers(T,M,h)

if isempty(T)
    T = M.T(~cellfun('isempty',M.T.F),:); 
end

for i=1:height(T)
    ph(i) = plot(M.t,T.X{i},'linew',2,'color',T.color(i,:)); hold on;
    ph(i).UserData = T.label{i};
end

set(gca,'YGrid','on','Box','off');

end