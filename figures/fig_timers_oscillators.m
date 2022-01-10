function [] = fig_timers_oscillators()

dbstop if error; close all;
h = fbmod_helpers;

fs = h.fs;
fs(end) = 20;

size_params = {...
    'radius',[1/4 1/10 1/6],...
    'fontsize',fs([2 3 2]),...
    'labeltype',{'full','full','full'}};

%%
M = fbmod_models('oscillators_ggg');

M.X0(M.typeixs('oscg')) = 0;
M = fbmod_add_action(M,{'extr_1','extr_1','extr_1'},{'oscg_1','oscg_2','oscg_3'},[0.010 0.015 0.020],1);
M = fbmod_add_action(M,{'etmr_2'},{'oscg_2'},M.tau(M.getid('etmr_2'),M.getid('gate_2')),-1);

M = fbmod_run_models(M);
M = fbmod_prep_schema(M,size_params{:});
M = h.adjpos(M,'osc_1',[0 0.5]);
M = h.adjpos(M,'osc_3',[0 0.5]);
M = h.adjpos(M,'extr_1',[1.20 1]);


%%
ax = stf([ones(5,2) repmat((2:6)',1,3)],...
    [0.01 .085 .01 .01],[0.15 0.01]);

%-------schema
axes(ax(1));
[H,C,h] = fbmod_draw_schema(M,h);

set(h.getconnh('extr_1_osc_1'),'linew',1);
set(h.getconnh('extr_1_osc_2'),'linew',1);
set(h.getconnh('extr_1_osc_3'),'linew',1);

connlabn = {'osc_1_osc_2','osc_2_osc_3','osc_3_osc_1'};
vas = {'bot','bot','top'};
has = {'lef','rig','cen'};
strs = h.latexify({'+','+','-'});
for i=1:length(connlabn)
    th(i) = label_connection(h.getconnh(connlabn{i}),[0.5 .5],strs{i},...
        'hori',has{i},'verti',vas{i},'fontsize',fs(end-1),'interpreter','latex');
end

%-------oscillator gates
axes(ax(2));
X = M.X(:,M.typeixs('oscg'));
for i=1:size(X,2)
    Gph(i) = plot(M.t,X(:,i),'-','linew',2,'color',M.G.color(i,:)); hold on;
end
set(gca,'YTick',[0 1]);

%-------oscillator amplitudes
axes(ax(3));
X = M.X(:,M.typeixs('oscr1'));
for i=1:size(X,2)
    Rph(i) = plot(M.t,X(:,i),'-','linew',2,'color',M.G.color(i,:)); hold on;
end
set(gca,'YTick',[0 1]);

%------oscillator activation
axes(ax(4));
Tosc = M.T(ismember(M.T.type,{'osc'}),:);
Oph = fbmod_plot_timers(Tosc,M,h);
set(gca,'YTick',[0 1]);

%------forces
axes(ax(5));
TX= M.T(ismember(M.T.type,{'osc' 'etmr' 'extr'}),:);
[Fph,Fth] = fbmod_plot_forces(TX,M,h); %#ok<*ASGLU>
set(Fth,'fontsize',fs(4));
for i=1:length(Fth)
    strs{i} = get(Fth(i),'string');
    set(Fth(i),'string',regexprep(strs{i},'\{\\rightarrow\}\{g_\d{1}\}',''));
end
ixs = find(contains(strs,'epsilon'));
set(Fth(ixs([1 end])),'visible','off');

%--------score
axes(ax(end));
h.score = fbmod_draw_score(M.G);


%%
set(ax,'box','off','tickdir','out','ticklen',0.005*[1 1],'fontsize',fs(end));
axrescaley(0.10,ax(2:end));
set(ax(2:end),'xlim',getlims(ax(2:end),'x'));
xlabel(ax(end),'time','fontsize',fs(2));
set(ax(2:end-1),'Xticklabel',[],'YGrid','on');

panlabs = {{'oscillator','gates'},{'oscillator','amplitudes'},{'oscillations'},'actions',{'activation','intervals'}};
stfig_panlab(ax(2:end),panlabs,'fontsize',fs(3),'xoff',-0.065,'yoff',0,'verti','top');

legend(Oph,Tosc.label,'fontsize',fs(3),'location','northeast','interpreter','latex','fontweight','normal');

%%
h.printfig(mfilename);


end

