function [] = fig_timers_actions_autonomy()

dbstop if error; close all;
h = fbmod_helpers;

fs = h.fs;
fs(3) = 20;

%% timer-gesture gate
M0 = fbmod_init_model({'g'},'dur',0.180);
M0 = fbmod_add_action(M0,{'etmr_1'},{'gate_1'},0.1,-1);

M0.tau(M0.tau==0.04) = 0.05;

M0 = fbmod_run_models(M0);

size_params = {'radius',[1/4 1/10 1/8],'fontsize',fs([2 3 2])};
M0 = fbmod_prep_schema(M0,size_params{:});

M0 = h.adjpos(M0,'etmr_1',[-0.1  0.25]);
M0 = h.adjpos(M0,'extr_1',[-0.75 -1]);
M0.T.label([2 4]) = {'${T}$','${\epsilon}$'};

%% autonomous vs. non-autonomous timers
M = fbmod_models('auto vs nonauto');
M = fbmod_run_models(M);

size_params = {'radius',[1/4 1/10 1/8],'fontsize',fs([2 3 2]),'labeltype',{'full' 'full' 'full'}};
M = fbmod_prep_schema(M,size_params{:});

M = h.adjpos(M,'extr_1',[0.875 0.25]);
M = h.adjpos(M,'extr_2',[-1.2 -0.5]);
M = h.adjpos(M,'etmr_2',[0  0.2]);

M.T.label = regexprep(M.T.label,'epsilon_\{1\}','epsilon^{\\prime}');
M.T.label = regexprep(M.T.label,'epsilon_\{2\}','epsilon_\{1\}');

%%

axpan0 = [ones(4,2) repmat((2:5)',1,3)];
axpan1 = [ones(3,2) repmat((2:4)',1,3)];
axpan = [axpan0; axpan1+max(axpan0(:))];
ax = stf(axpan,[0.05 .125 .01 .025],[0.20 0.01],'aspect',0.95);

ax0 = ax(1:5);
ax1 = ax(6:end);

for i=1:length(ax1)
    ax1(i).Position(2) = ax1(i).Position(2) - 0.05;
end

ax0(end).Position([2 4]) = ax0(end).Position([2 4]) + [0.05 -0.05];

%%

%-------schema
axes(ax0(1));
[H,C,h] = fbmod_draw_schema(M0,h);

%gate
rf = 3;
GATE = M0.G;
GATE.pos = GATE.pos + (GATE.radius)*[cos(3*pi/4) sin(3*pi/4)];
GATE.radius = GATE.radius/rf;
GATE.color = [1 .2 .2];
GATE.label = {'${G}$'};
GATE.facealpha = 1;
GATE.fontsize = fs(3);

H1 = h.getnode('extr_1');
H2 = h.getnode('etmr_1');

Gh.name = 'gate_1';
Gh.handles = draw_node(GATE);

C = [connection_constructor(H1,Gh); 
     connection_constructor(H2,Gh)]; 

%generate lines and scale
rl = [0.015 0.015]; %amount to reduce in normalized figure units
[rdux,rduy] = nfu2axu(rl(1),rl(2));

for i=1:height(C)
    cc = connecting_line(C.obj1(i),C.obj2(i));
    cc_len = sqrt(sum(diff(cc).^2));
    rd_len = sqrt(sum(rdux^2+rduy^2));
    sc = (cc_len-rd_len)/cc_len;
    ccsc = scale_connection(cc,sc);
    C.conn{i} = ccsc;
end

Ch = draw_connection(C,'linewidth',2);

delete(h.getconnh('extr_1_gest_1'));
delete(h.getconnh('etmr_1_gest_1'));

cprops = {'fontsize',fs(end),'interpreter','latex'};

Tix = @(name)find(ismember(M0.T.name,name));
 
str = h.spacer(h.latexify({sprintf('\\tau:%1.3f',M0.T.tau{Tix('extr_1')})}));
th(1) = label_connection(Ch{1},[0.5 0.25],str,cprops{:},'hori','right');

str = h.spacer(h.latexify({sprintf('\\tau_1:%1.3f',M0.T.tau{Tix('etmr_1')}(1))}));
th(2) = label_connection(Ch{2},[0.75 0.75],str,cprops{:},'hori','right');


%--------timer states
TX0 = M0.T(~cellfun('isempty',M0.T.F),:);
axes(ax0(2));
Th0 = fbmod_plot_timers(TX0,M0,h); %#ok<*NASGU>
set(Th0,'Linewidth',3);

%--------timer actions
axes(ax0(3));
[Ph,Fh] = fbmod_plot_forces(TX0,M0,h); %#ok<*ASGLU>
axrescaley(0.25);
set(Ph,'Linewidth',3);
set(Fh,'Fontsize',fs(3));

%---------gate states
axes(ax0(4));
gh = fbmod_plot_gates(M0,{},h);
gh.Color = GATE.color(1,:);
axrescaley(0.5);

%--------score
axes(ax0(5));
h.score = fbmod_draw_score(M0.G);

%------------------
ax0t = ax0(2:end);
set(ax0t,'xlim',getlims(ax0t,'x'));
set(ax0t(1:end-1),'xticklabel',[]);

panlabs = {{'timer','states'},{'actions'},{'gestural','gate'},{'activation','interval'}};
stfig_panlab(ax0t,panlabs,'fontsize',fs(3),'xoff',-0.30,'yoff',0,'verti','top','hori','left','fontweight','normal');

%%
%-------schema
axes(ax1(1));
[H,C,h] = fbmod_draw_schema(M,h);

lh = [h.getconnh('etmr_2_gest_2') h.getconnh('gest_2_etmr_2')];

harc1 = line2arcpair(lh,0.01*diff(xlim),0,0.02*diff(xlim),0,...
    'arrowprops',{'length',6,'tipangle',30},'linewidth',2);

set(gca,'Visible','off');

%--------timers
TX = M.T(~cellfun('isempty',M.T.F),:);
axes(ax1(2));
Th = fbmod_plot_timers(TX,M,h);

%---------actions
axes(ax1(3));
[ph,Fh] = fbmod_plot_forces(TX,M,h);
set(Fh,'Fontsize',fs(3));
axrescaley(0.35);

%--------score
axes(ax1(4));
h.score = fbmod_draw_score(M.G);

%----------
ax1t = ax1(2:end);
set(ax1t,'xlim',getlims(ax1t,'x'));
set(ax1t(1:end-1),'xticklabel',[]);

panlabs = {{'timer','states'},{'actions'},{'activation','intervals'}};
stfig_panlab(ax1t,panlabs,'fontsize',fs(3),'xoff',-0.30,'yoff',0,'verti','top','hori','left','fontweight','normal');


%%
set(ax,'box','off','tickdir','out','ticklen',0.005*[1 1],'fontsize',fs(end));

stfig_panlab([ax0(1) ax1(1)],{'A' 'B'},'fontsize',fs(2),'xoff',0,'yoff',0);

xlabel([ax1(end)],'time','fontsize',fs(end));

set(ax1(2),'Ytick',[0 max(M.tau(:))],'fontsize',fs(end)); 

strs = TX0.label(1:2);
legend(Th0,strs,'fontsize',fs(end),'location','northwest','interpreter','latex');

strs = TX.label(1:2);
strs{2} = '${\epsilon^{\prime}, \epsilon_1}$';
legend(Th(1:2),strs,'fontsize',fs(end),'location','northwest','interpreter','latex');


%%
h.printfig(mfilename);


end

