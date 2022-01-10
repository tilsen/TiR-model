function [] = fig_competitive_selection()

dbstop if error; close all;
h = fbmod_helpers;

%idea: extrinsic TiRs
N_sets = 3;
Ng_perset = 3;

%% ------- competitive queuing model of CV words

sn = @(str,n)sprintf('%s_%i',str,n);
sysnames = arrayfun(@(c){sprintf('C_%i',c),sprintf('V_%i',c),sprintf('R_%i',c)},(1:N_sets),'un',0);

utterance = {sysnames};

[M,U,P] = fbmod_gen_utterance_model('standard',utterance,'dur',1.25,'color',viridis(N_sets*Ng_perset)');
M = fbmod_run_models(M);

%%
fs = h.fs;
fs = [30 24 16 12];

schema_params = {...
    'radius',[1/3 1/8 1/5],...
    'fontsize',fs([2 4 3]),...
    'labeltype',{'full' 'minimal' 'full'}};

M = fbmod_prep_schema(M,schema_params{:});

M = h.adjpos(M,M.S.name(ismember(M.S.type,'extr')),[0 0.5]);
M.T.radius(ismember(M.T.type,'extr')) = 2*M.T.radius(ismember(M.T.type,'extr'));
M.G.row = (1+cell2mat(cellfun(@(c)mod((1:numel(c))-1,3),sysnames,'un',0)))';

%%
ax = stf([1 1 2 2 3 4 5]',[0.135 0.075 0.01 0.01],[0 0.02],'aspect',1.35);

conn_props = {...
    'reduce_lines',[0.005 0.005],...
    'linewidth',1,...
    'arrow_props',{'length',4,'tipangle',30}};

axes(ax(1));
H = fbmod_draw_schema(M,h,conn_props{:});

H = table2struct(H);

%rename mu systems
c=1;
for i=1:length(H)
    if regexp(H(i).name,'^extr_\d+$')
        if strcmp(H(i).handles.fh.Visible,'on')
            H(i).handles.th.String = ['${\mu_{' num2str(c) '}}$'];
            H(i).handles.th.FontSize = fs(2);
            c=c+1;
        end
    end
end


CQG = M.S(ismember(M.S.name,U.cqg_names),:);
CQ = M.S(ismember(M.S.name,U.cq_names),:);
CQ.lab = arrayfun(@(c){['\mu_' num2str(c)]},(1:height(CQ))');
colors = M.color';
cq_colors = colors(U.cq_gest_id,:);

%-------cq gates
axes(ax(end-2));
for j=1:height(CQG)
    y = M.X(:,CQG.id(j));
    fh(j) = fill([M.t; flipud(M.t)],[y; zeros(size(y))]-j,cq_colors(j,:),'facealpha',0.5,'edgecolor','none'); hold on;
end
axis tight;
set(gca,'YTick',[]);

%-------cq potential
axes(ax(end-3));

[D,E] = gen_cq_config(M.t,M.X(:,CQ.id),CQ.lab,'output_epochs','all');
E = format_cq_systems(E,'textloc','inside','colors',cq_colors,'fontsize',fs(3));
L = gen_cq_levels(D);
axcq = stfig_subaxpos(gca,[1:length(E)],[0 0 0 0 0 0]);

hCQ = draw_cq_potential(axcq,E,L,'steplinewidth',2);
ylims = getlims(axcq,'y');
set(axcq,'ylim',ylims);

%-------cq
axes(ax(end-1));
for j=1:height(CQ)
    y = M.X(:,CQ.id(j));
    phg(j) = plot(M.t,y,'color',cq_colors(j,:),'linew',2); hold on;
end
axis tight;
ylim([0 1.05]);
set(gca,'ygrid','on');

%------score
axes(ax(end));
h_sc = fbmod_draw_score(M.G);

%%
set(ax,'Box','off','tickdir','out','ticklen',0.001*[1 1],'fontsize',h.fs(end));
set(ax(2:end),'Xlim',minmax(M.t'),'xtick',[0:0.1:max(M.t)]);

set(ax(1:end-1),'XTickLabel',[]);
xlabel(ax(end),'time','fontsize',fs(2));

%stfig_panlab(ax,arrayfun(@(c){char(c+64)},(1:length(ax))),'verti','top','fontsize',h.fs(1),'xoff',-0.05);
panlabs = {'',...
    {'activation','potentials'},...
    {'\mu system','gates'},...
    {'\mu system','activation'},...
    {'gestural','activation','intervals'}};
stfig_panlab(ax,panlabs,'verti','top','fontsize',fs(2),'xoff',-0.05,'hori','right');

set(ax(2),'Visible','off');


%%
h.printfig(mfilename);

end
