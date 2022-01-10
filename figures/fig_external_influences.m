function [] = fig_external_influences()

dbstop if error; close all;
h = fbmod_helpers;
h.fs(end) = 14;

axpan = [repmat([1 1 1]',1,3) repmat([2 2 2]',1,4) repmat([3 4 5]',1,3)];
axpan = [axpan; nan(1,size(axpan,2)); axpan+max(axpan(:))];

ax = stf(axpan,[0.01 0.075 0.01 0.075],[0.1 0.01],'aspect',1.25);

ax1 = ax(1:5);
ax2 = ax(6:end);

control_transition(ax1,h);
frequency_variation(ax2,h);

%%
h.printfig(mfilename);


end

%%
function [] = frequency_variation(ax,h)

load([h.sim_dir 'fbmod_frequency_variation.mat']); %#ok<LOAD>

M{1} = h.adjpos(M{1},'extr_1',[1 0.5]);

lambda = R.lambda;
omegas = R.omega;

dg = R.deltas;
cc = lines(4);
colors = cc(3:4,:);

%%
ms = 24;

axes(ax(1));
[H,C,h] = fbmod_draw_schema(M{1},h);

axes(ax(2));
yyaxis left;
ph(1) = plot(lambda,dg(:,1)','k.-','linew',3,'markersize',ms); hold on;
ph(2) = plot(lambda,sum(dg,2),'.-','linew',3,'markersize',ms,'color',[.5 .5 .5]); hold on;

yyaxis right;
ph(3) = plot(lambda,omegas/(2*pi),'.-','linew',2,'color',colors(1,:),'markersize',ms); hold on;

Mix = [1 ceil(length(M)/2) length(M)];
for i=1:3
    axes(ax(i+2));
    fbmod_draw_score(M{Mix(i)}.G);
    text(max(M{1}.t),max(ylim)-0.1,sprintf('$\\lambda=%1.1f$',lambda(Mix(i))),...
        'hori','right','fontsize',h.fs(3),'verti','top','interp','latex');
end
axrescaley(0.025,ax(3:end));

%%
set(ax,'fontsize',h.fs(end),'tickdir','out','ticklen',0.005*[1 1]);
set(ax(3:end),'xlim',[0 M{1}.t(end)],'Xtick',0:0.1:max(M{1}.t));
set(ax([3 4]),'xticklabel',[]);
set(ax(2),'Box','off','XGrid','on','YGrid','on');

xlabel(ax(end),'time','fontsize',h.fs(3));

axes(ax(2));
yyaxis left;
ylabel('$\delta (\mathrm{s})$','fontsize',h.fs(2),'interp','latex');
set(gca,'YColor','k');

yyaxis right;
ylabel('$freq. (\mathrm{Hz})$','fontsize',h.fs(2),'interp','latex');
set(gca,'YColor','k');
axrescalex(0.05);
axrescaley(0.05);

xlabel('$\lambda$','fontsize',h.fs(2),'interp','latex');

legend(ph,{'$\delta_{CV}$','$\delta_{CR}$','$f$'},'fontsize',h.fs(3),'location','west','interpreter','latex');

stfig_panlab(ax(1),{'D'},'fontsize',h.fs(2),'xoff',0.1);
stfig_panlab(ax(2),{'E'},'fontsize',h.fs(2));
stfig_panlab(ax(3),{'F'},'fontsize',h.fs(2),'yoff',-0.1,'xoff',-0.01);
end

%% 
function [] = control_transition(ax,h)

load([h.sim_dir 'fbmod_control_transition.mat']);
M{1} = h.adjpos(M{1},'extr_1',[0 -0.5]);
lambda = R.lambda;
deltas = R.deltas;

TRix = find(R.stmr_first==1,1,'first');
cc = lines(4);
colors = cc(3:4,:);

%%

ms = 24;

axes(ax(1));
[H,C,h] = fbmod_draw_schema(M{1},h);

axes(ax(2));
yyaxis left;
ylim(minmax(deltas')+[-1 1]*0.025); hold on;
XX = [lambda(1) lambda(TRix)-diff(lambda(1:2))/2 lambda(end)];
fill(XX([1 1 2 2]),[ylim fliplr(ylim)],colors(1,:),'facealpha',0.2,'edgecolor','none'); hold on;
fill(XX([2 2 3 3]),[ylim fliplr(ylim)],colors(2,:),'facealpha',0.2,'edgecolor','none');
ph(1) = plot(lambda,deltas,'k.-','linew',3,'markersize',ms); hold on;

yyaxis right;
ph(2) = plot(R.lambda,R.alpha_etmr,'.-','linew',2,'color',colors(1,:),'markersize',ms); hold on;
ph(3) = plot(R.lambda,R.alpha_stmr,'.-','linew',2,'color',colors(2,:),'markersize',ms);

ixs = [1 find(R.stmr_first>0,1,'first') length(M)];
for i=1:length(M)
    axes(ax(i+2));
    fbmod_draw_score(M{i}.G);
    text(max(M{1}.t),max(ylim),sprintf('$\\lambda=%1.1f$',R.lambda(ixs(i))),...
        'hori','right','fontsize',h.fs(3),'verti','top','interp','latex');
end

%%
set(ax,'fontsize',h.fs(end),'tickdir','out','ticklen',0.005*[1 1]);
set(ax(3:end),'xlim',[0 M{1}.t(end)],'Xtick',0:0.1:max(M{1}.t));
set(ax([3 4]),'xticklabel',[]);
set(ax(2),'Box','off','XGrid','on','YGrid','on');

xlabel(ax(end),'time','fontsize',h.fs(3));

axes(ax(2));
yyaxis left;
ylabel('$\delta (s) = c_1 - V_1$','fontsize',h.fs(2),'interp','latex');
set(gca,'YColor','k');

yyaxis right;
ylabel('$\alpha$','fontsize',h.fs(2),'interp','latex');
set(gca,'YColor','k');
axrescalex(0.05);
axrescaley(0.05);

xlabel('$\lambda$','fontsize',h.fs(2),'interp','latex');

legh = legend(ph,{'$\delta$','$\alpha_{\hat{T_1}}$','$\alpha_{\bar{T_1}}$'},'fontsize',h.fs(3),'location','north','interpreter','latex');
legh.Position(2) = legh.Position(2) + 0.05;

stfig_panlab(ax(1),{'A'},'fontsize',h.fs(2),'xoff',0.1);
stfig_panlab(ax(2),{'B'},'fontsize',h.fs(2));
stfig_panlab(ax(3),{'C'},'fontsize',h.fs(2),'yoff',-0.1,'xoff',-0.01);


end

