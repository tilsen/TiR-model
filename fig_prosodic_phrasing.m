function [] = fig_prosodic_phrasing()

dbstop if error; close all;
h = fbmod_helpers;

fs = [30 24 16 12];

load([h.sim_dir 'fbmod_prosodic_phrasing.mat'],'M','U','P'); %#ok<*UNRCH>

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
Ux = {{{'\omega','Al'} {'\omega','and Bo'} {'\omega','or Cam'} {'\omega','were there'}};...
    {'iP' 'iP' 'iP'};
    {'IP' 'IP'};
    {'Utt'}};

% U = {{'\omega' '\omega' '\omega' '\omega'};...
%     {'iP' 'iP' 'iP'};
%     {'IP' 'IP'};
%     {'Utt'}};

S = [];
c=1;
for i=1:size(Ux,1)
    ux = Ux{i};
    for j=1:length(ux)
        S(c).id = c;
        S(c).lab = ux{j};
        S(c).row = i;
        S(c).col = j;
        S(c).xpos = nan;
        S(c).fontsize = 24;
        S(c).conn = nan;
        c=c+1;
    end
end

S = struct2table(S);

%positions
S.ypos = S.row;
S.xpos(1:4) = linspace(-2,2,4)';

%connections
S = repmat({S},1,2);


C{1} = cell2mat({[1 5],[2 5],[3 6],[4 7],[5 8],[6 8],[7 9],[8 10],[9 10]}');
C{2} = cell2mat({[1 5],[2 6],[3 6],[4 7],[5 8],[6 8],[7 9],[8 10],[9 10]}');

for j=1:length(C)
    for i=1:height(S{j})-1
        S{j}.conn(i) = C{j}(C{j}(:,1)==i,2);
    end
    
    %connection dependent positions
    for i=5:height(S{j})
        ixs = find(S{j}.conn==i);
        S{j}.xpos(i) = mean(S{j}.xpos(ixs));
    end
end

%%
ax = stf([repmat([1 2],2,1) [3; 4]],[0.02 0.02 0.02 0.02],[0.05 0.10],'aspect',3);

for i=1:2
    axes(ax(i));
    s = S{i};
    for j=1:height(s)
        
        if ~isnan(s.conn(j))
            xx = [s.xpos(j) s.xpos(s.conn(j))];
            yy = [s.ypos(j) s.ypos(s.conn(j))];
            plot(xx,yy,'k-','linew',2); hold on;
        end
        
        th(i,j) = text(s.xpos(j),s.ypos(j),s.lab{j},...
            'hori','center','verti','mid','fontsize',s.fontsize(j),'BackGroundcolor','w',...
            'margin',0.01,'Fontname','cambria');
    end
    axis tight;
    
    p1 = th(i,2).Extent*[1 0 1 0; 0 1 0 0]' - [-0.01 0.1];
    p0 = p1-[0 .25];
    arrow(p0,p1,'linewidth',3,'tipangle',30,'length',6,'color','r');
    
end

lims = getlims(ax(1:2),'xy');
set(ax(1:2),'xlim',lims(1,:),'ylim',lims(2,:));
axrescalex(0.1,ax(1:2));
axrescaley([-0.1 0.2],ax(1:2));
set(ax(1:2),'Visible','off');


%----------------
axcq1 = stfig_subaxpos(ax(3),1:height(CQ1),[0 0 0 0 0 0]);
axcq0 = stfig_subaxpos(ax(4),1:height(CQ1),[0 0.1 0 0.1 0 0]);

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
    
    axcq0s{i} = stfig_subaxpos(axcq0(i),1:height(CQ0),[0.001 0 0.001 0 0 0]);
    [D0,E0] = gen_cq_config(M.t,M.X(:,CQ0.id),CQ0.lab,'output_epochs','anyselected');
    E0 = format_cq_systems(E0,'textloc','above','colors',cq0_colors,'fontsize',fs(3));
    L0 = gen_cq_levels(D0);
    
    set(axcq0(i),'XTick',[],'YTick',[],'Box','on');
    
    hCQ0 = draw_cq_potential(axcq0s{i},E0,L0,'steplinewidth',2);
    for j=1:length(hCQ0)
        set(hCQ0(j).sys_th(ishandle(hCQ0(j).sys_th)),'Visible','off');
    end
end

ax0s = [axcq0s{:}];
set(ax0s,'Visible','off');
axrescalex(0.25,ax0s(1));

%
axbak = stbgax;
P0 = [vertcat(axcq1.Position)*[1 0 0.5 0; 0 1 0 0.20]']; 
P1 = [vertcat(axcq0.Position)*[1 0 0.5 0; 0 1 0 1]'];
for i=1:length(axcq1)
    arh = arrow(P0(i,:),P1(i,:),'length',8,'tipangle',30);
    th = label_connection(arh,0.5,sprintf('$\\lambda_%i$',i),'fontsize',36,'interpreter','latex','xoff',0);
    
    %if i==2, th.Color = 'r'; th.FontWeight = 'bold'; end
end

%%
set(ax,'Visible','off');
panlabs = {'A','B','C','D'};
stfig_panlab(ax,panlabs,'verti','top','fontsize',fs(1),'xoff',-0.01,'hori','right');

%%
h.printfig(mfilename);

end
