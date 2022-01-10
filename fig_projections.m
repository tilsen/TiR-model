function []  = fig_projections()

dbstop if error; close all;
h = fbmod_helpers;

%% prepare models

%CV from babbling
M = fbmod_models({'hybrid sensory','hybrid internal'});
M = fbmod_run_models(M);

M{1}.color(:,2) = mean(M{1}.color(:,[1 3]),2);
M{2}.color = viridis(size(M{2}.color,2))';

fs = [36 28 18 16];
size_params = {...
    'radius',[1/3 1/4 1/4],...
    'fontsize',fs([3 4 4]),...
    'labeltype',{'full' 'minimal' 'full'}};

M = fbmod_prep_schema(M,size_params{:});

M{1}.G.row = [1 2 3 1 3]';
M{2}.G.row = [1 2 3 1 3]';

%x: developmental time
%y: utterance time
%z: selection scale

for i=1:length(M)
    H(i).M = M{i};
    H(i).M.T.has_action(ismember(H(i).M.T.type,'extr')) = false;
end

H(1).pos = [1 1 1];
H(2).pos = [2 1 2];

H(1).spos = [1.5 M{1}.G.t_on{1} 1; ...
             1.5 M{2}.G.t_on{4} 1];
H(2).spos = [2.5 M{2}.G.t_on{2}+0.050 2];

H(1).scolor = [M{1}.color(:,1)'; ...
               M{1}.color(:,4)'];
H(2).scolor =  M{2}.color(:,3)';

%%
ax = stf([1 1],[0.05 0.05 0.01 0.01]);

%----- setup axes
view(3);
xlim([0.5 3]); ylim([0 0.800]); zlim([0 3]); hold on;
set(gca,'XGrid','on','YGrid','on','ZGrid','on','GridAlpha',0.25,...
    'YDir','reverse','TickDir','in','TickLen',0.001*[1 1],'Box','on');
light('position',[0 1 3],'style','local');

%[az,el] = view;
rots = [11 -18 90];

xlabh = xlabel('developmental time \rightarrow','rotation',rots(1),'fontsize',fs(1));
ylabh = ylabel('utterance time \rightarrow','rotation',rots(2),'fontsize',fs(1));
zlabh = zlabel('internalization \rightarrow','rotation',rots(3),'fontsize',fs(1));

xlabh.Position = xlabh.Position + [0.65    0   0.05];
ylabh.Position = ylabh.Position + [0 -.25   0.025];

set(gca,'XTickLabel',[],'ZTickLabel',[],'XTick',[1 2 3],'ZTick',[1 2]);

%----gestural scores
for i=1:length(H)
    H(i).score = draw_score(H(i));    
end

%----selection system spheres
for i=1:length(H)
    H(i).spheres = draw_sphere(H(i));    
end

%----projection
for i=1:length(H)
    H(i).projections = draw_projections(H(i));
    H(i).projlines = draw_projection_lines(H(i));
end

%----schemas
for i=1:length(M)
    
    H(i).schema = fbmod_draw_schema(H(i).M,h,'format_axes',false,'linewidth',1);
    reposition_schema(H(i));   
end

%----utterance scope
for i=1:length(M)
    H(i).uttscope = draw_utterance_scope(H(i)); 
end

%----structure:
o11 = H(1).projections(1);
o12 = H(1).projections(2);
o2 = H(2).projections;

H(1).structline(1) = line([o11.XData o2.XData],[o11.YData o2.YData],[o11.ZData o2.ZData],'color','k','linew',2);
H(1).structline(2) = line([o12.XData o2.XData],[o12.YData o2.YData],[o12.ZData o2.ZData],'color','k','linew',2);

text(o11.XData,o11.YData,o11.ZData,'\mu','fontsize',fs(2),'hori','center','verti','top');
text(o12.XData,o12.YData,o12.ZData,'\mu','fontsize',fs(2),'hori','center','verti','top');
text(o2.XData,o2.YData,o2.ZData,'\sigma','fontsize',fs(2),'hori','center','verti','bot');

%%

h.printfig(mfilename);

end


%%
function [] = write_animated_gif(filename,FR,delaytime)

for idx = 1:length(FR)
    im = frame2im(FR(idx));
    [A,map] = rgb2ind(im,256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',delaytime);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',delaytime);
    end
end
end

%% draw score
function [sc] = draw_score(H)
G = H.M.G;
pos = H.pos;
xh = 0.75;
xpos = linspace(pos(1),pos(1)+xh,numel(G.row)+1);
for i=1:height(G)
    xx = xpos(G.row(i)+[0 0 1 1]);
    yy = [G.t_on{i} G.t_off{i}];
    yy = yy([1 2 2 1]);
    zz = zeros(size(xx));
    sc(i) = patch(xx,yy,zz,G.color(i,:),'facealpha',G.facealpha(i),'Edgecolor','k'); hold on;
end
end

%% draw sphere
function [sph] = draw_sphere(H)

pos = H.spos;

r = [0.055 0.02 0.125];
[xx,yy,zz] = sphere(200);
xx = r(1)*xx;
yy = r(2)*yy;
zz = r(3)*zz;

for i=1:size(pos,1)
    C = H.scolor(i,:);
    sph(i) = surf(xx+pos(i,1),yy+pos(i,2),zz+pos(i,3),'FaceColor',C,...
        'EdgeColor','none','SpecularStrength',0.15);
end
end

%% draw projections
function [proj] = draw_projections(H)

sph = H.spheres;
for i=1:length(sph)
    y = mean(sph(i).YData(:));
    z = mean(sph(i).ZData(:));
    proj(i) = plot3(max(zlim),y,z,'k.','markersize',30);
end

end

%% reposition schemas on x-z plane
function [] = reposition_schema(H)

schema = H.schema;
pos = H.pos;
%reposition and rescale x-coordinates

%x-offset of schema
hh = [schema.hvec{:}];
hh_patch = hh(ismember(get(hh,'Type'),'patch'));
xlims = cell2mat(arrayfun(@(c){minmax(c.XData(:)')},hh_patch'));
ylims = cell2mat(arrayfun(@(c){minmax(c.YData(:)')},hh_patch'));
xlims = minmax(xlims(:)');
ylims = minmax(ylims(:)');
xrsc = 1/(1.5*diff(xlims));
yrsc = 1/(diff(ylims));

for i=1:height(schema)
    h = schema.hvec{i};
    for j=1:length(h)
        switch(h(j).Type)
            case 'patch'
                h(j).XData = h(j).XData - xlims(1);
                h(j).XData = h(j).XData*xrsc;
                h(j).XData = h(j).XData + pos(1);
                h(j).YData = h(j).YData - ylims(1);
                h(j).YData = h(j).YData*yrsc;
                h(j).YData = h(j).YData + pos(3);                
            case 'text'
                h(j).Position(1) = h(j).Position(1) - xlims(1);
                h(j).Position(1) = h(j).Position(1)*xrsc;
                h(j).Position(1) = h(j).Position(1) + pos(1);
                h(j).Position(2) = h(j).Position(2) - ylims(1);
                h(j).Position(2) = h(j).Position(2)*yrsc;
                h(j).Position(2) = h(j).Position(2) + pos(3);                
        end
    end
end

for i=1:height(schema)
    switch(schema.type{i})
        case 'gesture'
            fh = schema.handles{i}.fh;
            th = schema.handles{i}.th;
            fh.ZData = fh.YData;
            fh.YData = min(ylim)*ones(size(fh.YData));
            th.Position = th.Position([1 3 2]);    
            
        case 'timer'
            fh = schema.handles{i}.fh;
            th = schema.handles{i}.th;
            fh.ZData = fh.YData;
            fh.YData = min(ylim)*ones(size(fh.YData));
            th.Position = th.Position([1 3 2]);             
            
        case 'sens'
            fh = schema.handles{i}.fh;
            th = schema.handles{i}.th;
            fh.ZData = fh.YData;
            fh.YData = min(ylim)*ones(size(fh.YData));
            th.Position = th.Position([1 3 2]);              
            
        case 'connection'
            arh = schema.hvec{i};
            arh.ZData = arh.YData;
            arh.YData = min(ylim)*ones(size(arh.YData));
            
    end
end

end

%% draw utterance scope
function [uttscope] = draw_utterance_scope(H)

sph = H.spheres;
G = H.M.G; 
G.t_on = cell2mat(G.t_on);
G.t_off = cell2mat(G.t_off);

if numel(sph)>1
    sphy = arrayfun(@(c)mean(c.YData(:)),sph)';
    dGon = pdist2(sphy,G.t_on);
    [~,ixs] = min(dGon);
else
    ixs = ones(1,height(G));
end

xrng = minmax(reshape([H.score.XData],1,[]));

for i=1:length(sph)
   
    %scope
    yy = [min(G.t_on(ixs==i)) max(G.t_off(ixs==i))];
    yy = yy([1 1 2 2 1]);
    xx = xrng([1 2 2 1 1]);
    
    %source
    ys = mean(sph(i).YData(:));
    xs = mean(sph(i).XData(:));
    zs = min(sph(i).ZData(:));
    
    for j=1:4
        uttscope(i,j) = patch([xs xx(j+(0:1)) xs],[ys yy(j+(0:1)) ys],[zs 0 0 zs],H.scolor(i,:),'facealpha',0.1);
    end
    
end


end

%% draw projection lines
function [lineh] = draw_projection_lines(H)

sph = H.spheres;
for i=1:length(sph)
    
    xs = min(sph(i).XData(:));
    ys = mean(sph(i).YData(:));
    zs = mean(sph(i).ZData(:));
    
    xx = H.projections(i).XData-0.01*diff(xlim);
    yy = H.projections(i).YData;
    zz = H.projections(i).ZData;
    
    %lineh(i) = line([xs xx],[ys yy],[zs zz],'color','k','linew',2);
    lineh(i) = arrow([xs ys zs],[xx yy zz],'color','k','linew',2);
    
end

end

%%
% 
% 
% dtoff = linspace(0,1,5); dtoff = dtoff(2:end-1); %centers of gestural systems in developmental time
% utoff = {dtoff,[mean(dtoff(1:2)) dtoff(3)],0.5}; %...utterance time
% scoff = [dtoff]; %...scale
% 
% TV = [1 0 0; 0 1 0; 0 0 1];
% %score times
% cw = 0.150;
% vw = 0.400;
% times{1} = [utoff{1}(1)+[-cw cw]/2; nan nan; nan nan];
% times{1}(2,:) = [times{1}(1,2)+[0 vw]];
% times{1}(3,:) = [times{1}(2,2)+[0 cw]];
% times{2} = [times{1}(1,:); times{1}(2,:)-0.100; nan nan];
% times{2}(3,:) = times{2}(2,2)+[0 cw];
% times{3} = [times{2}(1:2,:); times{2}(3,:)-0.100];
% 
% STR{1} = {'X','X','X'};
% STR{2} = {'\mu','\mu'};
% STR{3} = {'\sigma'};
% 
% cmap = jet(10);
% colors{1} = cmap(2:4,:);
% colors{2} = cmap(6:7,:);
% colors{3} = cmap(9,:);
% 
% fsl = 26;
% labrots = [12.5, -20];
% toz= 0.025;
% tox = -0.005;
% rots = labrots;
% fs = 26;
% zfac = 2.05; %stretches spheres
% textzboost = 0.015;
% 
% %%
% axes(ax(1));
% 
% Y{1} = [dtoff(1) utoff{1}(1) scoff(1); dtoff(1) utoff{1}(2) scoff(1); dtoff(1) utoff{1}(3) scoff(1)]';
% Y{2} = [dtoff(2) utoff{2}(1) scoff(2); dtoff(2) utoff{2}(2) scoff(2)]';
% Y{3} = [dtoff(3) utoff{3}(1) scoff(3)]';
% [xx,yy,zz] = sphere(200);
% 
% for j=1:length(Y)
%     X = Y{j};
%     A = [1 0 0; 0 0 0; 0 0 1];
%     B = [0 0 0; 0 1 0; 0 0 1];
%     C = [1 0 0; 0 1 0; 0 0 1];
%     Xa = A*X;
%     Xb = B*X;
%     Xb(1,:) = 1-Xb(1,:);    
%     XA{j} = Xa;
%     XB{j} = Xb;
% end
% 
% gray = [.5 .5 .5];
% Cg = squeeze(permute(repmat(0.85*ones(1,3),[1 1 size(zz,1) size(zz,2)]),[3 4 2 1]));
% 
% %draw surfaces
% cc = 0.025;
% for j=1:length(Y)
%     X = Y{j};
%     for i=1:size(X,2)
%         C = squeeze(permute(repmat(colors{j}(i,:),[1 1 size(zz,1) size(zz,2)]),[3 4 2 1]));
%         sh{j,i} = surf(X(1,i)+cc*xx,X(2,i)+cc*yy,X(3,i)+zfac*cc*zz,C,'EdgeColor','none','SpecularStrength',0.15); hold on;
%         sh_gray{j,i} = surf(X(1,i)+cc*xx,X(2,i)+cc*yy,X(3,i)+zfac*cc*zz,Cg,'EdgeColor','none','SpecularStrength',0.15,'visible','off'); hold on;
% 
%     end
% end
% 
% % draw projected structures
% hs_xoff = 0;
% for i=1:2
%     switch(i)
%         case 1
%             hsprops = {'color','k','linew',1};
%             
%         case 2
%             hsprops = {'color','k','linew',2,'visible','off'};
%     end
%     hs(1,i) = line([1 1]+hs_xoff,[XB{3}(2,1) XB{2}(2,1)],[XB{3}(3,1) XB{2}(3,1)],hsprops{:});
%     hs(2,i) = line([1 1]+hs_xoff,[XB{3}(2,1) XB{2}(2,2)],[XB{3}(3,1) XB{2}(3,2)],hsprops{:});
%     hs(3,i) = line([1 1]+hs_xoff,[XB{2}(2,1) XB{1}(2,1)],[XB{2}(3,1) XB{1}(3,1)],hsprops{:});
%     hs(4,i) = line([1 1]+hs_xoff,[XB{2}(2,1) XB{1}(2,2)],[XB{2}(3,1) XB{1}(3,2)],hsprops{:});
%     hs(5,i) = line([1 1]+hs_xoff,[XB{2}(2,2) XB{1}(2,3)],[XB{2}(3,2) XB{1}(3,3)],hsprops{:});
% end
% 
% %draw projection
% for j=1:length(Y)
%     X = Y{j}; Xa = XA{j}; Xb = XB{j};
%     strs = STR{j};
%     
%     for i=1:size(X,2)
%         uprojph(j,i) = plot3([X(1,i) Xa(1,i)],[X(2,i) Xa(2,i)],[X(3) Xa(3,i)],'Color','k','linestyle',':','linew',2);
%         uprojlh(j,i) = plot3(Xa(1,i),Xa(2,i),Xa(3,i),'ko','markerfaceColor','k'); hold on;
%         
%         dprojph(j,i) = plot3([X(1,i) Xb(1,i)],[X(2,i) Xb(2,i)],[X(3,i) Xb(3,i)],'Color','k','linestyle',':','linew',2);
%         dprojlh(j,i) = plot3(Xb(1,i),Xb(2,i),Xb(3,i),'ko','markerfaceColor','k'); hold on;
%     end
% end
% set(ax,'XGrid','on','YGrid','on','ZGrid','on','YDir','reverse');
% xlim([0 1]);
% ylim([0 1]);
% zlim([0 1]);
% ticks = linspace(0,1,32);
% zticks = linspace(0,1,16);
% 
% set(gca,'XGrid','on','YGrid','on','ZGrid','on','YDir','reverse',...
%     'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],...
%     'XTick',ticks,'YTick',ticks,'ZTick',zticks,'TickLen',[0.01 0.01]);
% 
% xlabel('developmental time \rightarrow','rotation',labrots(1));
% ylabel('utterance time \rightarrow','rotation',labrots(2));
% zlabel('selection scale \rightarrow','rotation',90);
% set(ax.XLabel,'verti','top','hori','center','Position',[0.5 1.05],'fontsize',fsl);
% set(ax.YLabel,'verti','top','hori','center','Position',[-0.05 0.5],'fontsize',fsl);
% set(ax.ZLabel,'verti','bottom','hori','center','Position',[-0.025 0 0.5],'fontsize',fsl);
% 
% 
% 
% 
% %draw surfaces text
% for j=1:length(Y)
%     X = Y{j}; strs = STR{j};
%     for i=1:size(X,2)
%          sth(j,i) = text(...
%             min(sh{j,i}.XData(:)),...
%             max(sh{j,i}.YData(:)),...
%             max(sh{j,i}.ZData(:))+textzboost,strs{i},'verti','bottom','hori','center','fontsize',fs); hold on;
%     end
% end
% 
% %draw projected structure text
% for j=1:length(Y)
%     X = Y{j}; Xa = XA{j}; Xb = XB{j}; strs = STR{j};
%     for i=1:size(X,2)
%         uprojth(j,i) = text(Xa(1,i),Xa(2,i),Xa(3,i)+toz,strs{i},'verti','bottom','hori','center','fontsize',fs,'rotation',rots(1)); hold on;
%         dprojth(j,i) = text(Xb(1,i)+tox,Xb(2,i),Xb(3,i)+toz,strs{i},'verti','bottom','hori','center','fontsize',fs,'rotation',rots(2)); hold on;
%     end
% end
% dprojth(1,3).Position = dprojth(1,3).Position + [-.01 .0225 0]; %fix for line overlap
% 
% w = 0.2;
% yr = [.05 .95];
% [xx,yy] = meshgrid(linspace(0-w/2,0+w/2,100),linspace(yr(1),yr(2),100));
% xx = xx+0.01*randn(size(xx));
% yy = yy+0.01*randn(size(yy));
% 
% %fields
% for i=1:3
%     nlh(i) = line(dtoff(i)+xx(:),yy(:),-.001*ones(1,length(xx(:))),'marker','.','linestyle','none','color',[.5 .5 .5],'markersize',2);
%     nfb(i) = line(dtoff(i)+[w w -w -w w]/2,yr([1 2 2 1 1]),-0.001*ones(1,5),'color',.75*ones(1,3));
% end
% set(gca,'CLipping','off');
% 
% %gestures
% gw = 0.05;
% for i=1:3
%     for j=1:size(TV,1)
%         gsfh(i,j) = fill3((dtoff(i)+(j-1.5)*gw)+[-gw -gw 0 0],[times{i}(j,[1 2 2 1])],zeros(1,4),[.5 .5 .5],'facealpha',0.2);
%     end
% end
% 
% %gesture projections
% xtimes{1} = times{1};
% xtimes{2} = [times{2}(1,1) times{2}(2,2); times{2}(3,:)];
% xtimes{3} = [times{3}(1,1) times{3}(3,2)];
% for i=1:3
%     si = sh(i,:);
%     si = si(cellfun(@(c)~isempty(c),si(:)));
%     for j=1:length(si)
%         ssxx = dtoff(i)+[-w -w w w]/2;
%         ssyy = [xtimes{i}(j,[1 2 2 1])];
%         ssfh(i,j) = fill3(ssxx,ssyy,[0 0 0 0]-0.0001,colors{i}(j,:),...
%             'facealpha',0.45,'edgecolor',[.4 .4 .4]);
%         ssfh_gray(i,j) = fill3(ssxx,ssyy,[0 0 0 0]-0.0001,gray,...
%             'facealpha',0.45,'edgecolor',[.4 .4 .4],'visible','off');        
%         
%         for k=1:length(ssxx)
%             gplh(i,j,k) = line([mean(sh{i,j}.XData(:)) ssxx(k)],[mean(sh{i,j}.YData(:)) ssyy(k)],[min(sh{i,j}.ZData(:)) 0],'color',colors{i}(j,:));
%             gplh_gray(i,j,k) = line([mean(sh{i,j}.XData(:)) ssxx(k)],[mean(sh{i,j}.YData(:)) ssyy(k)],[min(sh{i,j}.ZData(:)) 0],'color',gray,'visible','off');
%            
%         end
%         
%         gpph1(i,j) = patch([mean(sh{i,j}.XData(:)) ssxx(1:2) mean(sh{i,j}.XData(:))],...
%             [mean(sh{i,j}.YData(:)) ssyy(1:2) mean(sh{i,j}.YData(:))],...
%             [min(sh{i,j}.ZData(:)) 0 0 min(sh{i,j}.ZData(:))],colors{i}(j,:),'edgecolor','none','facealpha',0.2);
%         
%         gpph1_gray(i,j) = patch([mean(sh{i,j}.XData(:)) ssxx(1:2) mean(sh{i,j}.XData(:))],...
%             [mean(sh{i,j}.YData(:)) ssyy(1:2) mean(sh{i,j}.YData(:))],...
%             [min(sh{i,j}.ZData(:)) 0 0 min(sh{i,j}.ZData(:))],gray,'edgecolor','none','facealpha',0.2,'visible','off');        
%         
%         gpph2(i,j) = patch([mean(sh{i,j}.XData(:)) ssxx(2:3) mean(sh{i,j}.XData(:))],...
%             [mean(sh{i,j}.YData(:)) ssyy(2:3) mean(sh{i,j}.YData(:))],...
%             [min(sh{i,j}.ZData(:)) 0 0 min(sh{i,j}.ZData(:))],colors{i}(j,:),'edgecolor','none','facealpha',0.2);
%         
%          gpph2_gray(i,j) = patch([mean(sh{i,j}.XData(:)) ssxx(2:3) mean(sh{i,j}.XData(:))],...
%             [mean(sh{i,j}.YData(:)) ssyy(2:3) mean(sh{i,j}.YData(:))],...
%             [min(sh{i,j}.ZData(:)) 0 0 min(sh{i,j}.ZData(:))],gray,'edgecolor','none','facealpha',0.2,'visible','off');       
%     end
% end
% 
% % ssfh = ssfh';
% % ssfh_gray = ssfh_gray';
% % gplh = permute(gplh,[2 1 3]);
% % gplh_gray = permute(gplh_gray,[2 1 3]);
% % gpph1 = gpph1';
% % gpph2 = gpph2';
% % gpph1_gray = gpph1_gray';
% % gpph2_gray = gpph2_gray';


%end

% update_objects('toggle_neurons');        %toggle selection set fields (expansion of scale)
% update_objects({'toggle_uprojections','toggle_uprojection_text','toggle_uprojection_dots'});   %toggle projection over developmental time
% update_objects({'toggle_gfields','toggle_sfields'});        %toggle selection set fields (expansion of scale)
% update_objects('toggle_labels');        %toggle selection set fields (expansion of scale)
% update_objects('toggle_axes');        %toggle selection set fields (expansion of scale)
% 
% print('-dpng','-r400','fig_projection_01.png');
% 
% update_objects('toggle_darkstructure');        %toggle selection set fields (expansion of scale)
% update_objects('toggle_axes');        %toggle selection set fields (expansion of scale)
% update_objects('toggle_dprojections');   %toggle projection over utterance time
% update_objects({'toggle_systems','toggle_system_text'});
% 
% print('-dpng','-r400','fig_projection_02.png');
% 
% update_objects('toggle_darkstructure');        %toggle selection set fields (expansion of scale)
% update_objects({'toggle_structure','toggle_structure_dots','toggle_structure_text'});
% update_objects({'toggle_systems','toggle_gfields','toggle_sfields'});
% update_objects({'toggle_axes','toggle_labels'});        %toggle selection set fields (expansion of scale)
% 
% print('-dpng','-r400','fig_projection_03.png');
% 
% update_objects({'toggle_structure','toggle_structure_dots','toggle_structure_text'});
% update_objects({'toggle_gfields','toggle_sfields'});
% update_objects({'toggle_axes'});        %toggle selection set fields (expansion of scale)
% update_objects('toggle_dprojections');
% 
% print('-dpng','-r400','fig_projection_04.png');
% 
% update_objects('CV_labels');
% update_objects('toggle_dprojections');
% update_objects({'toggle_uprojection_text','toggle_uprojection_dots'});
% update_objects({'toggle_gfields','toggle_sfields'});
% 
% set(gcf,'Color','w');
% update_objects('focus_1'); FR(1) = getframe(gcf);
% update_objects('focus_2'); FR(2) = getframe(gcf);
% update_objects('focus_3'); FR(3) = getframe(gcf);
% 
% write_animated_gif('fig_organization_anim.gif',FR,1.5);
% 
% update_objects('focus_2'); 
% print('-dpng','-r400','fig_projection_05.png');
% 
% update_objects('focus_3'); 
% print('-dpng','-r400','fig_projection_06.png');
% 
% clear('FR');
% update_objects('focus_1'); FR(1) = getframe(gcf);
% update_objects('focus_dprojections_1'); FR(end+1) = getframe(gcf);
% update_objects('clear_dprojections');
% update_objects('focus_2'); FR(end+1) = getframe(gcf);
% update_objects('focus_dprojections_2'); FR(end+1) = getframe(gcf);
% update_objects('clear_dprojections');
% update_objects('focus_3'); FR(end+1) = getframe(gcf);
% update_objects('focus_dprojections_3'); FR(end+1) = getframe(gcf);
% 
% write_animated_gif('fig_projection_anim.gif',FR,0.75);
% 
%     function [] = update_objects(option,varargin)
%         
%         G = false(3);
%         offon = {'off','on','Visible'};
%         Xcolorknone = {'k','none','Xcolor'};
%         Ycolorknone = {'k','none','Ycolor'};
%         Zcolorknone = {'k','none','Zcolor'};
%         toggleprop = @(c,d)d{mod(strcmp(get(c,d{3}),d{2})+1,2)+1};
%         gethandles = @(c)c(arrayfun(@(d)ishandle(d),reshape(c(:),1,[])));
%         gethandlesc = @(c)c(cellfun(@(c)~isempty(c),reshape(c,1,[])));
%         
%         if iscell(option)
%             for ii=1:length(option)
%                 update_objects(option{ii})
%             end
%             return;
%         end
%         
%         switch(option)
%             
%             case 'toggle_axes'
%                 set(ax,'XColor',toggleprop(ax,Xcolorknone));
%                 set(ax,'YColor',toggleprop(ax,Ycolorknone));
%                 set(ax,'ZColor',toggleprop(ax,Zcolorknone));
%                 
%             case 'toggle_labels'          
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),[ax.XLabel ax.YLabel ax.ZLabel]);
%                 set([ax.XLabel ax.YLabel ax.ZLabel],'color','k');
%                 
%             case 'toggle_neurons'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),[nlh nfb]);
%                 
%             case 'toggle_sfields'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles([ssfh(:); gplh(:); gpph1(:); gpph2(:)]));
%                 
%             case 'toggle_gfields'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(gsfh(:)));
% 
%             case 'toggle_dprojections'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(dprojph(:)));
%                 
%             case 'focus_dprojections_1'
%                 G = logical([1 1 1; 0 0 0; 0 0 0]);
%                 arrayfun(@(c)set(c,'Visible','on'),gethandles(dprojph(G(:)))); 
%                 arrayfun(@(c)set(c,'Visible','off'),gethandles(dprojph(~G(:))));
%                 
%             case 'focus_dprojections_2'
%                 G = logical([0 0 0; 1 1 0; 0 0 0]);
%                 arrayfun(@(c)set(c,'Visible','on'),gethandles(dprojph(G(:)))); 
%                 arrayfun(@(c)set(c,'Visible','off'),gethandles(dprojph(~G(:))));  
%                 
%             case 'focus_dprojections_3'
%                 G = logical([0 0 0; 0 0 0; 1 0 0]);
%                 arrayfun(@(c)set(c,'Visible','on'),gethandles(dprojph(G(:)))); 
%                 arrayfun(@(c)set(c,'Visible','off'),gethandles(dprojph(~G(:))));      
%                 
%             case 'clear_dprojections'
%                 arrayfun(@(c)set(c,'Visible','off'),gethandles(dprojph(:)));                    
%                 
%             case 'toggle_uprojections'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(uprojph(:)));
%                 
%             case 'toggle_uprojection_dots'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(uprojlh(:)));
%                 
%             case 'toggle_uprojection_text'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(uprojth(:)));
%                 
%             case 'toggle_structure_dots'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(dprojlh(:)));
% 
%             case 'toggle_structure_text'
%                  arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(dprojth(:)));
% 
%             case 'toggle_structure'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(hs(:,1)));
%                 
%             case 'toggle_darkstructure'      
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(hs(:,2)));
%                 
%             case 'toggle_system_text'
%                 arrayfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandles(sth(:)));
%                 
%             case 'toggle_systems'
%                 cellfun(@(c)set(c,'Visible',toggleprop(c,offon)),gethandlesc(sh));
%                 
%             case 'CV_labels'
%                 dlabs = {'C' 'V' 'C'; '\mu' '\mu' ''; '\sigma' '' ''}; dlabs = gethandlesc(dlabs);
%                 ulabs = {'\{C\}\{V\}\{C\}' '' ''; '\{CV\}\{C\}' '' ''; '\{CVC\}' '' ''}; 
%                 arrayfun(@(c,d)set(d,'String',dlabs{c}),1:numel(gethandles(dprojth)),gethandles(dprojth));
%                 arrayfun(@(c,d)set(d,'String',ulabs{c}),1:numel(gethandles(uprojth)),gethandles(uprojth));
% 
%             case 'X_labels'
%                 dlabs = {'X' 'X' 'X'; '\mu' '\mu' ''; '\sigma' '' ''};
%                 ulabs = {'\{X\}' '\{X\}' '\{X\}'; '\{\mu\}' '\{\mu\}' ''; '\{\sigma\}' '' ''};
%                 arrayfun(@(c,d)set(d,'String',dlabs{c}),1:numel(gethandles(dprojth)),gethandles(dprojth));
%                 arrayfun(@(c,d)set(d,'String',ulabs{c}),1:numel(gethandles(uprojth)),gethandles(uprojth));
%                     
%             case 'focus_3'
%                 go = [1 1; 1 2; 1 3; 2 1; 2 2]; G(go(:,1),go(:,2)) = true;
%   
%             case 'focus_2'
%                 go = [1 1; 1 2; 1 3; 3 1]; G(go(:,1),go(:,2)) = true;             
%                 
%             case 'focus_1'
%                 go = [2 1; 2 2; 3 1]; G(go(:,1),go(:,2)) = true;
%                 
%             case 'focus_0'
%                 G = false(3);
%          end
%         
%         switch(option)
%             case {'focus_0','focus_1','focus_2','focus_3'}
%                arrayfun(@(c)set(c,'Visible','off'),gethandles([gpph1_gray(:); gpph2_gray(:); ssfh_gray(:); gpph1(:); gpph2(:); ssfh(:); gplh(:); gplh_gray(:)]));
%                 cellfun(@(c)set(c,'Visible','off'),gethandlesc([sh(:) sh_gray(:)]));
%                 arrayfun(@(c)set(c,'Visible','on'),gethandles([gpph1_gray(G(:)); gpph2_gray(G(:)); ssfh_gray(G(:)); gplh_gray(reshape(repmat(G(:),1,1,4),[],1))]));
%                 cellfun(@(c)set(c,'Visible','on'),gethandlesc(sh_gray(G(:))));
%                 arrayfun(@(c)set(c,'Visible','on'),gethandles([gpph1(~G(:)); gpph2(~G(:)); ssfh(~G(:)); gplh(reshape(repmat(~G(:),1,1,4),[],1))]));
%                 cellfun(@(c)set(c,'Visible','on'),gethandlesc(sh(~G(:))));
%                 arrayfun(@(c)set(c,'fontweight','normal','color',[.5 .5 .5]),gethandles([uprojth(G(:)); dprojth(G(:))]));
%                 arrayfun(@(c)set(c,'fontweight','bold','color',[0 0 0]),gethandles([uprojth(~G(:)); dprojth(~G(:))]));
%                 arrayfun(@(c)set(c,'markerface',[.5 .5 .5],'markersize',6),gethandles([uprojlh(G(:)); dprojlh(G(:))]));
%                 arrayfun(@(c)set(c,'markerface',[0 0 0],'markersize',8),gethandles([uprojlh(~G(:)); dprojlh(~G(:))]));                
%         end
%        
%     end
% end


