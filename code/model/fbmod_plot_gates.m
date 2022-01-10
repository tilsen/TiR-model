function [ph] = fbmod_plot_gates(M,names,h)

if isempty(names)
    names = M.S.name(ismember(M.S.type,'gate'));
end

for i=1:length(names)
    Six = find(ismember(M.S.name,names{i}));
    X(:,i) = M.X(:,Six);
    color(i,:) = M.color(:,M.S.type_id(Six))';
end

for i=1:length(names)
    ph(i) = plot(M.t,X(:,i),'linew',2,'color',color(i,:)); hold on;
end

set(gca,'YTick',[0 1],'YTickLabel',{'clo' 'op'},'YGrid','on');
axis tight;
axrescaley(0.05);

end