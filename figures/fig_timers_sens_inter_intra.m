function [] = fig_timers_sens_inter_intra()

dbstop if error; close all;
h = fbmod_helpers;

fs = h.fs;
fs(3) = 20;

size_params = {...
    'radius',[1/4 1/10 1/6],...
    'fontsize',fs([2 3 2]),...
    'labeltype',{'full' 'full' 'full'}};

%% external vs internal 
M0 = fbmod_models('external vs internal');
M0 = fbmod_run_models(M0);
M0 = fbmod_prep_schema(M0,size_params{:});
M0 = h.adjpos(M0,'extr_1',[0 -0.5]);

%% timer-gesture gate
M1 = fbmod_models('inter vs intra');
M1 = fbmod_run_models(M1);
M1 = fbmod_prep_schema(M1,size_params{:});
M1 = h.adjpos(M1,'extr_1',[0 -0.5]);

%%
axpan = [ones(3,2) repmat((2:4)',1,3)];
axpan = [axpan; axpan+max(axpan(:))];
ax = stf(axpan,[0.05 .075 .02 .05],[0.10 0.01],'aspect',0.95);

ax0 = ax(1:4);
ax1 = ax(5:end);

for i=1:length(ax1)
    ax1(i).Position(2) = ax1(i).Position(2) - 0.05;
end

ax0(end).Position([2 4]) = ax0(end).Position([2 4]) + [0.05 -0.05];
ax1(end).Position([2 4]) = ax1(end).Position([2 4]) + [0.05 -0.05];

%%

%-------schema
axes(ax0(1));
[H,C,h] = fbmod_draw_schema(M0,h);

%--------timer states
axes(ax0(2));
Th0 = fbmod_plot_timers({},M0,h); %#ok<*NASGU>

%--------timer actions
axes(ax0(3));
[Ph,Fh] = fbmod_plot_forces({},M0,h); %#ok<*ASGLU>
axrescaley(0.25);
set(Fh,'Fontsize',fs(3));

%--------score
axes(ax0(4));
h.score = fbmod_draw_score(M0.G);

%------------------
ax0t = ax0(2:end);
set(ax0t,'xlim',getlims(ax0t,'x'));
set(ax0t(1:end-1),'xticklabel',[]);

% panlabs = {{'timer','states'},{'actions'},{'gestural','gate'},{'activation','interval'}};
% stfig_panlab(ax0t,panlabs,'fontsize',fs(3),'xoff',-0.30,'yoff',0,'verti','top','hori','left','fontweight','normal');

%%
%-------schema
axes(ax1(1));
[H,C,h] = fbmod_draw_schema(M1,h);

%--------timers
axes(ax1(2));
Th1 = fbmod_plot_timers({},M1,h);

%---------actions
axes(ax1(3));
[ph,Fh] = fbmod_plot_forces({},M1,h);
set(Fh,'Fontsize',fs(3));
axrescaley(0.35);

%--------score
axes(ax1(4));
h.score = fbmod_draw_score(M1.G);

%----------
ax1t = ax1(2:end);
set(ax1t,'xlim',getlims(ax1t,'x'));
set(ax1t(1:end-1),'xticklabel',[]);

% panlabs = {{'timer','states'},{'actions'},{'activation','intervals'}};
% stfig_panlab(ax1t,panlabs,'fontsize',fs(3),'xoff',-0.30,'yoff',0,'verti','top','hori','left','fontweight','normal');


%%
set(ax,'box','off','tickdir','out','ticklen',0.005*[1 1],'fontsize',fs(end));

stfig_panlab([ax0(1) ax1(1)],{'A' 'B'},'fontsize',fs(2),'xoff',0,'yoff',0);

xlabel([ax1(end)],'time','fontsize',fs(end));

set(ax0t(1),'Ytick',[0 max(M0.tau(:))],'fontsize',fs(end)); 
set(ax1t(1),'Ytick',[0 max(M1.tau(:))],'fontsize',fs(end)); 

legend(Th0,{Th0.UserData},'fontsize',fs(end),'location','northwest','interpreter','latex');

legend(Th1,{Th1.UserData},'fontsize',fs(end),'location','northwest','interpreter','latex','NumColumns',2);

%%
h.printfig(mfilename);


end

