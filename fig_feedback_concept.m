function [] = fig_feedback_concept()

dbstop if error; close all;
h = fbmod_helpers();

%% define and run model
M = fbmod_init_model({'g_1','g_2'},'dur',0.500,'freq',5.0,'dt',0.0005);

M = fbmod_add_action(M,{'extr_1','extr_2'},{'gate_1','gate_2'},[0 0.001],1); 
M = fbmod_add_action(M,'stmr_1','gate_1',0.250,-1);
M = fbmod_add_action(M,'etmr_2','gate_2',0.248,-1);

M = fbmod_run_models(M);

%% format schema
fs = [48 36 28 18];
size_params = {'radius',[1/3 1/8 1/3],'fontsize',fs([1 3 1])};

M = fbmod_prep_schema(M,size_params{:});

M.Y.pos(:,2) = M.Y.pos(:,2)-0.5;
M.Y.label{1} = regexprep(M.Y.label{1},'g','g^{\\prime}');

% exclude extrinsic timers:
M.T.has_action(ismember(M.T.type,'extr')) = false;

%% plot
ax = stf([1 1],[0.01 0.01 0.01 0.01],'aspect',1.0);

%-------schema
[H,C,h] = fbmod_draw_schema(M,h);

lh = [h.getconnh('etmr_2_gest_2') h.getconnh('gest_2_etmr_2')];

harc1 = line2arcpair(lh,0.01*diff(xlim),0,0.02*diff(xlim),0,...
    'arrowprops',{'length',6,'tipangle',30},'linewidth',2);

g1 = h.getnodeh('gest_1'); g1 = g1(1);
s1 = h.getnodeh('sens_1'); s1 = s1(1);

PP = [min(g1.XData)-0.03 mean(g1.YData);...
    nan(2,2); ...
    min(s1.XData)-0.03 mean(s1.YData)];

PP(2,:) = mean(PP([1 end],:)) + [-0.25 0.55];
PP(3,:) = mean(PP([1 end],:)) + [-0.25 -0.55];

harc2 = draw_arc(PP,'linewidth',2);

yo = 0.15;
har1 = h.getconnh('sens_1_stmr_1');

label_connection(harc1(1).lh,[0.01 0.5],{'internal','feedback'},'fontsize',fs(2),'verti','mid','xoff',0.15);
label_connection(har1,[0.5 0.5],{'external','feedback'},'fontsize',fs(2),'verti','bot','hori','left','xoff',0.05);

xx = [-.5 .5];
PP = [xlim+xx fliplr(xlim+xx); (max(s1.YData)+yo)*[1 1] (max(ylim)+0.1)*[1 1]]';
PP(end+1,:) = PP(1,:);
ph = plot(PP(:,1),PP(:,2),'k--','linewidth',2);

yo = 0.05;
text(max(xlim)-0.05,ph.YData(1)+yo,'central nervous system','hori','right','verti','bot','fontsize',fs(3));
text(max(xlim)-0.05,ph.YData(1)-yo,{'peripheral nervous system / ','environment'},'hori','right','verti','top','fontsize',fs(3));

%%
set(ax,'visible','off')


%%
h.printfig(mfilename);



end
