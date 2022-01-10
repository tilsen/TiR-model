function [] = fig_threebody_predictions()

dbstop if error; close all;
h = fbmod_helpers;

fs = h.fs;
cortype = 'cor';

%% load models and run (just for schemas)
models = {
    'internal_ggg_1'
    'oscillators_ggg'
    'external_ggg'
    'internal_ggg_2'
    'internal_ggg_0'
    'internal_ggg_3'}';

model_descr = {...
    'shared trigger',...
    'coupled oscillators',...
    'external feedback chain',...
    'internal feedback chain',...
    'independent ext. triggers',...
    'hybrid independent/shared'};

M = fbmod_models(models);
M = fbmod_run_models(M);

%% prepare schemas
size_params = {...
    'radius',[1/4 1/10 1/6],...
    'fontsize',fs([3 4 3]),...
    'labeltype',{'full' 'short' 'full'}};

M = fbmod_prep_schema(M,size_params{:});

%oscillator trigger
osc_mod = find(ismember(models,'oscillators_ggg'));
oth_mod = setdiff(1:numel(models),osc_mod);

M{osc_mod} = h.adjpos(M{osc_mod},'extr_1',[1 0]);

for i=oth_mod
    M{i} = h.adjpos(M{i},'extr_1',[0 -0.5]);
    M{i} = h.adjpos(M{i},'extr_2',[0 -0.5]);
    M{i} = h.adjpos(M{i},'extr_3',[0 -0.5]);
end

M{3} = h.adjpos(M{3},'extr_1',[0 -0.25]);

%% load simulation results
load([h.sim_dir 'fbmod_threebody_covariance.mat'],'R');

R.noise_global(R.parset==1) = R.omega_glob(R.parset==1);
R.noise_local(R.parset==1) = R.omega_local(R.parset==1);
R.noise_global(R.parset==2) = R.alpha_glob(R.parset==2);
R.noise_local(R.parset==2) = R.alpha_local(R.parset==2);

for i=1:3
    R.N_on(:,i) = cellfun(@(c)numel(c),R.t_on(:,i));
end

%exclusions:

R = R(~any(R.N_on~=1,2),:);

R.t_on = cell2mat(R.t_on);
R.don = diff(R.t_on,[],2);

%% some "errors" must be excluded:
R.don_neg = R.don<0;

% crosstab(R.model,R.don_neg(:,1));
% crosstab(R.model,R.don_neg(:,2));

R = R(all(R.don>0,2),:);

%%

for i=1:length(models)
    ix = ismember(R.model,models{i});
    noise_global = unique(R.noise_global(ix));
    noise_local = unique(R.noise_local(ix));
    R.noise_global_ix(ix) = arrayfun(@(c)find(ismember(noise_global,c)),R.noise_global(ix));
    R.noise_local_ix(ix) = arrayfun(@(c)find(ismember(noise_local,c)),R.noise_local(ix));
end

Rs = grpstats(R,{'model' 'noise_global_ix' 'noise_local_ix'},{'@(c){c}' 'mean' 'std'},'DataVars',{'don'});

Rs.corr = cellfun(@(c,d)corr(c,d),Rs.Fun1_don(:,1),Rs.Fun1_don(:,2));

Rs = sortrows(Rs,{'model' 'noise_global_ix' 'noise_local_ix'});

Np = numel(unique(R.noise_global_ix));

for i=1:length(models)
    c = Rs.corr(ismember(Rs.model,models{i}));
    C(i,:,:) = reshape(c,Np,[])';
end

%%
Nm = length(models);
axpan = reshape(1:6,2,[])';
ax = stf([axpan Nm+axpan],...
    [0.05 0.075 0.01 0.05],[0.05 0.05]);

for i=(Nm+1):2:length(ax)
    ax(i).Position(1) = ax(i).Position(1) + 0.035;
end
for i=1:Nm
    ax(i).Position(2) = ax(i).Position(2) - 0.02;
end

ms = {'s','o','o','o','o'};
ls = {'-','-','-','-','-'};
lw = [2 2 2 3 3];

colors = stf_colormap(Np,[0 0 0],[1 0 0]);

C(:,1,1) = nan;

arprops = {'length',4,'tipangle',30};

for i=1:Nm
    axes(ax(i));
    [H,Cn,h] = fbmod_draw_schema(M{i},h,'arrow_props',arprops,'linewidth',1);
    
    set(gcf,'currentaxes',ax(i+Nm));
    for a=1:size(C,2)
        ph(i,a) = plot(squeeze(C(i,a,:)),ms{a},'color',colors(a,:),'linestyle',ls{a},'linew',lw(a),'markerfacecolor','w'); hold on;
    end
    
end

%%
axC = ax(Nm+1:end);
axis(axC,'tight');
set(axC,'Xlim',[0.75 Np+0.25]);

switch(cortype)
    case 'cov'
        set(axC,'ylim',getlims(axC,'y'));
    case 'cor'
        set(axC,'ylim',[-1 1]);
end

noise_labs = [{'0'} repmat({''},1,Np-2) {'max'}];
noise_labsf = [{'0'} {'+'} {'++'} {'+++'} {'max'}];

axrescaley(0.10,axC);
set(ax,'Box','off','YGrid','on','fontsize',h.fs(end),'tickdir','out');
set(axC,'XTick',[1:Np],'XTickLabel',noise_labs);

xlabel(ax([end-1 end]),'local noise strength');

set(axC(2:2:end),'YTickLabel',[])
set(axC(1:end-2),'XTickLabel',[])

arrayfun(@(c)plot3(get(c,'xlim'),[0 0],[-1 -1],'--','color',[.5 .5 .5],'linew',1,'parent',c),axC);
arrayfun(@(c)plot3(get(c,'xlim'),[-1 -1],[-1 -1],'-','color',[.5 .5 .5],'linew',2,'parent',c),axC);
arrayfun(@(c)plot3(get(c,'xlim'),[1 1],[-1 -1],'-','color',[.5 .5 .5],'linew',2,'parent',c),axC);

legh = legend(ph(1,:),noise_labsf,...
    'fontsize',h.fs(end),'interpreter','none','location','southwest','NumColumns',2);

legh.Title.String = 'global noise';

pl = arrayfun(@(c){char(c+64)},(1:Nm)');
panlabs = cellfun(@(c,d){['{\bf{' c ')}} ' d]},pl,model_descr');

ax(Nm-1).Position(2)= ax(Nm-1).Position(2)-0.05;
ax(Nm).Position(2)= ax(Nm).Position(2)-0.05;

stfig_panlab(ax(1:Nm),panlabs,'fontsize',fs(3),'hori','left','fontweight','normal',...
    'xoff',0,'yoff',[-0.01 -0.01 -0.04 -0.04 -0.04 -0.04],'verti','bot','positionmode','outerposition');

stfig_panlab(ax(Nm+1:end),cellfun(@(c){[c '^{\prime}']},pl),'fontsize',fs(3),'hori','left','verti','bot','yoff',-0.025,'xoff',-0.01);

ax(Nm-3).Position(2)= ax(Nm-3).Position(2)-0.025;
ax(Nm-2).Position(2)= ax(Nm-2).Position(2)-0.025;

%%
h.printfig(mfilename);

end
