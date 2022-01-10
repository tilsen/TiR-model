function [] = fbmod_quickplot(M)

dbstop if error;
h = fbmod_helpers();
fs = h.fs;

if ~isfield(M,'G')
    M = fbmod_prep_schema(M);
end

%types of systems to plot
types = unique(M.S.type);
types = setdiff(types,'oscr2');

types_gates = {'cqg','exg','gate','oscg'};

Npanels = length(types)+2+1; %add two panels for forces and one for score

%% system dynamics
Nc = 4;
ax = stf(reshape(1:16,[],Nc),...
    [0.05 0.05 0.01 0.01],[0.05 0.01],'handlearray','matrix');

t = M.t;
colors = M.color';

%--------system states
for i=1:length(types)
    set(gcf,'currentaxes',ax(i));
    set(gca,'xlim',minmax(M.t'),'XGrid','on'); hold on;   
    ixs = ismember(M.S.type,types{i});
    X = M.X(:,ixs);
    names = M.S.name(ixs);
    
    switch(types{i})
        case types_gates
            
            imagesc(t,(1:length(names)),X'); hold on; colormap([1 1 1; 1 .5 .5]); 
            axis tight; set(gca,'YDir','reverse'); drawnow;
            glh = stfig_gridlines('y','linew',0.5,'color',[.5 .5 .5]);
            for j=1:length(names)
                t_first = t(find(X(:,j)>0,1,'first'));
                if isempty(t_first), t_first = t(1); end                
                th{i}(j) = text(t_first,j,names{j},'verti','mid','fontsize',12,'interpreter','none','color','k');
            end                         
        otherwise
            for j=1:size(X,2)
                ph{i}(j) = plot(t,X(:,j),'linew',2,'color',colors(j,:)); hold on;
            end
            legend(ph{i},names,'location','southeast','NumColumns',2,'interpreter','none');
    end
    ylabel(types{i},'fontsize',fs(3));
    
    ax(i).UserData = types{i};
end

%--------forces
TT = {M.T(ismember(M.T.type,{'osc' 'extr'}),:),...
      M.T(~ismember(M.T.type,{'osc' 'extr'}),:)};

for j=1:length(TT)
    set(gcf,'currentaxes',ax(end-3+j));
    [phF{j},thF{j}] = fbmod_plot_forces(TT{j},M,h); %#ok<*NASGU,*AGROW>
    set(gca,'UserData','');
end

%-------score
axes(ax(end));
M.G.fontsize = 12*ones(height(M.G),1);
H = fbmod_draw_score(M.G);
set(gca,'UserData','');

%%
axis(ax(1:end-3),'tight');
set(ax,'xlim',minmax(M.t'),'XGrid','on');

axrs = ax(~ismember({ax.UserData},types_gates));
axrescaley(0.05,axrs);

set(ax,'ticklen',0.003*[1 1],'tickdir','out','box','off');

axm = reshape(ax,[],Nc);
set(axm(1:end-1,:),'XTickLabel',[]);

set(ax(1:end-1),'YGrid','on');

end