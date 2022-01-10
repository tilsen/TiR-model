function [] = fig_multilevel_competitive_selection()

dbstop if error; close all;
h = fbmod_helpers;

fs = [30 24 16 12];

load([h.sim_dir 'fbmod_multilevel_competitive_selection.mat'],'M','U','P'); %#ok<*UNRCH>

%%
pwrds = unique(U.pwrd);
for i=1:length(pwrds)
    g = U.gests(U.pwrd==i);
    g = [g{:}];
    pwrd_num_gests(i) = numel(g);
end

rows = arrayfun(@(c){(1:c)'},pwrd_num_gests)';
M.G.row = vertcat(rows{:});

%CQG1 = M.S(ismember(M.S.name,P.cqg_name),:);
CQ1 = M.S(ismember(M.S.name,P.cq_name),:);
CQ1.lab = arrayfun(@(c){['\mathcal{C}_' num2str(c)]},(1:height(CQ1))');
colors = M.color';
cq1_colors = colors(P.last_gest_typeid,:);

%%
ax = stf([1 2 3]',[0.135 0.075 0.01 0.01],[0 0.02],'aspect',1.35);

%------score
axes(ax(end));
h_sc = fbmod_draw_score(M.G);
for i=1:length(h_sc)
    h_sc(i).th.Position(1) = h_sc(i).th.Position(1)-0.02;
end

axcq1 = stfig_subaxpos(ax(1),1:height(CQ1),[0 0 0 0 0 0]);
axcq0 = stfig_subaxpos(ax(2),1:height(CQ1),[0 0 0 0 0.01 0]);


%-------cq level 1 potential
[D1,E1] = gen_cq_config(M.t,M.X(:,CQ1.id),CQ1.lab,'output_epochs','anyselected');
E1 = format_cq_systems(E1,'textloc','inside','colors',cq1_colors,'fontsize',fs(3));
L1 = gen_cq_levels(D1);
hCQ1 = draw_cq_potential(axcq1,E1,L1,'steplinewidth',2);
ylims = getlims(axcq1,'y');
set(axcq1,'ylim',ylims);

%-------cq level 0 potentials
for i=1:height(CQ1)
    
    %CQG0 = M.S(ismember(M.S.name,U.cqg_names(U.pwrd==i)),:);
    CQ0 = M.S(ismember(M.S.name,U.cq_names(U.pwrd==i)),:);
    
    nums = U.name(U.pwrd==i);
    nums = cellfun(@(c)c(end),nums);
    
    CQ0.lab = arrayfun(@(c){['\mu_{' c '}']},nums);
    cq0_colors = colors(U.cq_gest_id(U.pwrd==i),:);
    
    axcq0s{i} = stfig_subaxpos(axcq0(i),1:height(CQ0),[0.01 0 0.01 0 0 0]);
    [D0,E0] = gen_cq_config(M.t,M.X(:,CQ0.id),CQ0.lab,'output_epochs','anyselected');
    E0 = format_cq_systems(E0,'textloc','above','colors',cq0_colors,'fontsize',fs(3));
    L0 = gen_cq_levels(D0);
    
    set(axcq0(i),'XTick',[],'YTick',[],'Box','on');
   
    hCQ0 = draw_cq_potential(axcq0s{i},E0,L0,'steplinewidth',2);
end

ax0s = [axcq0s{:}];
set(ax0s,'Visible','off');

axrescalex(0.25,axcq0s{3});

%%
set(ax,'Box','off','tickdir','out','ticklen',0.001*[1 1],'fontsize',h.fs(end));
set(ax(1:end-1),'Visible','off');
set(ax(end),'Xlim',minmax(M.t'),'xtick',[0:0.25:max(M.t)]);

xlabel(ax(end),'time (s)','fontsize',fs(2));

panlabs = {...
    {'concept','system','competitive','selection'};
    {'\mu system','competitive','selection'};
    {'gestural','activation','intervals'}};

stfig_panlab(ax,panlabs,'verti','top','fontsize',fs(2),'xoff',-0.01,'hori','right');

%%
h.printfig(mfilename);

end
