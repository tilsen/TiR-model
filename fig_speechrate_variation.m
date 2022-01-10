function [] = fig_speechrate_variation()

dbstop if error; close all;
h = fbmod_helpers;

load([h.sim_dir 'fbmod_speechrate_boundary_variation.mat'],'R','M0');

fs = h.fs;
fs(2) = 24;
fs(3) = 16;

size_params = {...
    'radius',[1/4 1/8 1/6],...
    'fontsize',fs([2 3 2]),...
    'labeltype',{'full' 'short' 'full'}};

M0 = fbmod_prep_schema(M0,size_params{:});

M0.T.label = regexprep(M0.T.label,'\epsilon','\mu');

extr_ixs = find(ismember(M0.S.type,'extr'));
for i=extr_ixs'
    M0 = h.adjpos(M0,M0.S.name{i},[1 0.5]);
end

M0.sets = cellfun(@(c)str2double(c),regexp(M0.name,'\d+','match'));

Nm = height(R);

%%

ix_word_ons = arrayfun(@(c)find(M0.sets==c,1,'first'),unique(M0.sets));
ix_word_off = arrayfun(@(c)find(M0.sets==c,1,'last'),unique(M0.sets));

ix_last_CR = find(ismember(M0.name,{'C_3' 'R_3'}));
ix_last_Rc = find(ismember(M0.name,{'R_3' 'c_3'}));
ix_last_cr = find(ismember(M0.name,{'c_3' 'r_3'}));

for i=1:Nm
    t_on = cell2mat(R.t_on(i,:));
    R.T_ON(i,1:length(t_on)) = t_on;
    
    t_off = cell2mat(R.t_off(i,:));
    R.T_OFF(i,1:length(t_off)) = t_off;    
end

R.word_dur = R.T_OFF(:,ix_word_off)-R.T_ON(:,ix_word_ons);
R.final_CR = diff(R.T_ON(:,ix_last_CR),[],2);
R.final_V  = diff(R.T_ON(:,ix_last_Rc),[],2);
R.final_cr = diff(R.T_ON(:,ix_last_cr),[],2);

R.freq = R.omega/(2*pi);

%%
axpan = [1 1 2; 3 4 5];
ax = stf(axpan,[0.05 0.085 0.05 0.01],[0.085 0.05]);

%ax(1).Position(3) = ax(1).Position(3)-0.02;

lprops = {'linew',2,'markerfacecolor','w'};

arprops = {'length',4};

%----schema
axes(ax(1));
[H,C,h] = fbmod_draw_schema(M0,h,...
    'arrow_props',arprops,'reduce_lines',[0.005 0.005], 'linewidth',1); 
 
%----params
axes(ax(2));
pcolors = [0.5 .5 .5; 0.5 0 0; 0 0 0];

yyaxis left;
parh(1) = plot(R.lambda,R.alpha_etmr,'o-','color',pcolors(1,:),lprops{:}); hold on;
            plot(R.lambda,R.alpha_etmr_last,'o:','color',pcolors(1,:),lprops{:}); hold on;
parh(2) = plot(R.lambda,R.alpha_stmr,'s-','color',pcolors(2,:),lprops{:}); 
            plot(R.lambda,R.alpha_stmr_last,'s:','color',pcolors(2,:),lprops{:}); hold on;
set(gca,'Ycolor','k');
ylabh(1) = ylabel('\alpha','fontsize',h.fs(2));
axis tight;

yyaxis right;
parh(3) = plot(R.lambda,R.freq,'o-','color',pcolors(3,:),lprops{:}); hold on;
set(gca,'Ycolor','k');
ylabh(2) = ylabel('freq. (Hz)','fontsize',h.fs(2));
axis tight;


%----onset times
axes(ax(3));
colors = M0.color';
for i=1:size(R.T_ON,2)
    ph(i) = plot(R.lambda,R.T_ON(:,i),'o-','color',colors(i,:),lprops{:}); hold on;
end

%----word durations
axes(ax(4));
wcolors = colors(ix_word_ons,:);
for i=1:size(R.word_dur,2)
    wh(i) = plot(R.lambda,R.word_dur(:,i),'o-','color',wcolors(i,:),lprops{:}); hold on;
end

%----CR vs. V vs. cr duration
axes(ax(5));
ff = {'final_CR','final_V','final_cr'};
icolors = lines(3);
for i=1:length(ff)
    ih(i) = plot(R.lambda,R.(ff{i}),'o-','color',icolors(i,:),lprops{:}); hold on;
end


%%
set(ax,'Fontsize',h.fs(end),'Box','off','tickdir','out','ticklen',0.003*[1 1]);

axes(ax(2));
yyaxis left; axis tight; axrescaley(0.05); 
yyaxis right; axis tight; axrescaley(0.05);
set(gca,'XTickLabel',[]);

axis(ax(3:end),'tight');
axrescalex(0.05,ax(2:end));
axrescaley(0.05,ax(3:end));
set(ax(2:end),'YGrid','on','XGrid','on');

set(ylabh,'fontsize',fs(2));

xlabel(ax(3:end),'\lambda','fontsize',h.fs(3));
ylabel(ax(3),'gestural onset time (s)','fontsize',h.fs(3));

ylabel(ax(3),'gestural onset time (s)','fontsize',h.fs(3));
ylabel(ax(4),'word dur. (s)','fontsize',h.fs(3));
ylabel(ax(5),'interval dur. (s)','fontsize',h.fs(3));

legend(parh,{'internal TiRs','sensory TiRs','osc. freq.'},'fontsize',fs(2),'location','northeast');
legend(wh,{'word 1','word 2','word 3'},'fontsize',fs(2),'location','northwest');
legend(ih,{'C3-R3','R3-c3','c3-r3'},'fontsize',fs(2),'location','northwest');

xoffs = -0.11*ones(size(ax));
xoffs(1) = 0;
stfig_panlab(ax,arrayfun(@(c){char(c+64)},(1:length(ax))),'verti','top','xoff',xoffs);


%%

h.printfig(mfilename);


end