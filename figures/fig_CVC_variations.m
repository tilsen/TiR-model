function [] = fig_CVC_variations()

dbstop if error; close all;
h = fbmod_helpers;

%% load and run models
mods = [{'maximally sensory'};{'maximally internal'};{'oscillator triggered'};{'hybrid model'}];
M = fbmod_models(mods);
MM = fbmod_run_models(M);

%% prep schemas
fs = [28 18 16 12];
size_params = {'radius',[1/3 1/9 1/3],'fontsize',fs([2 3 2])};

for i=1:length(MM)
    M = MM{i};   
    M = fbmod_prep_schema(M,size_params{:});
   
    %system layout adjustments
    switch(M.model_descr{:})
        case {'oscillator triggered','hybrid model'}
%             M = adjpos(M,'osc_2',[0.25 0]);            
%             M = adjpos(M,'osc_1',[0.5 0]);    
%             M = adjpos(M,'osc_3',[0 0]);    
%             M = adjpos(M,'osc_4',[0 0.25]);    
%             M = adjpos(M,'osc_5',[0 0.25]);      
            M.T.has_action(ismember(M.T.type,'extr')) = false;
    end
       
    M.T.label(ismember(M.T.type,'itmr')) = regexprep(M.T.label(ismember(M.T.type,'itmr')),'{T}_','');
    
    MM{i} = M;
    
end

%% plot
Nr = length(MM);
ax = stf([2 2],[0.01 .01 .01 .025],[0.05 0.05]);

axpos = cell2mat(arrayfun(@(c){get(c,'position')},ax'));

xw = 0.65;
yw = 0.75;
for i=1:size(axpos,1)
    axsc(i) = axes('position',...
        [axpos(i,1) + xw*axpos(i,3) axpos(i,2) + yw*axpos(i,4) (1-xw)*axpos(i,3)   (1-yw)*axpos(i,4)]);
end

ylims = [-1.75 3.25];

for i=1:Nr
    axes(ax(i));
    h.H{i} = fbmod_draw_schema(MM{i},h);
    axis tight; axis 'equal'; 
    ylim(ylims);
        
    axes(axsc(i));
    MM{i}.G.row = [1 2 3 1 3]';
    h.score{i} = fbmod_draw_score(MM{i}.G);
    set(gca,'Visible','off');
end

%%
labs = cellfun(@(c){c.model_descr{:}},MM);
plabs = arrayfun(@(c){char(c+64)},(1:length(labs))');
labs = cellfun(@(c,d){[c ') ' d]},plabs,labs');

for i=1:length(ax)
    ax(i).Position(1)=ax(i).Position(1)-0.05;
end

stfig_panlab(ax,labs,'xoff',-0.1,'fontsize',fs(1),'yoff',-0.05,'hori','left');

%%
h.printfig(mfilename);

end

