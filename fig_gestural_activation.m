function [] = fig_gestural_activation()

dbstop if error; close all;
h = fbmod_helpers;
load([h.sim_dir 'CVC_example_trajectories.mat']);

basepan = [1 2; 1 2; 1 3; 1 4];
axpan = [repmat(basepan(:,1),1,4) repmat(basepan(:,2),1,4)];
ax = stf(axpan,[0.025 0.085 0.01 0.065],[0.075 0.075],'aspect',2.0);
fs = [36 28 20 18]; h.fs = fs;

%%
set(gcf,'currentaxes',ax(1));
N = base_system(h);
[~,N] = draw_node(N);
axis tight; xlim(xlim+[-1 0]);
axis equal; hold on; 

C=[];
for i=1:height(N)
    conn = N.conn{i};
    for j=1:length(conn)
        C = [C; node_connector(N(i,:),N(conn(j),:),'reduceline',[0.015 0.015])];
    end
end

%treat segments differently
C.linestyle = repmat({'-'},height(C),1);
C.linestyle(end-2:end) = {'--'};
C.conntype(end-2:end) = {'none'};

C = draw_connections(C);

text(min(xlim),3,{'gestural','systems'},'fontsize',fs(2));
text(min(xlim),2,{'vocal tract','systems'},'fontsize',fs(2));
text(min(xlim),1,{'articulator','systems'},'fontsize',fs(2));

N.fh(end).Visible = 'off';
N.fh(end-1).Visible = 'off';

%%
arprops = {'length',20,'tipangle',30,'baseangle',90};
heights = 0.20;
tvlabs = {'LA','PHAR'};

xlims = [-.2 .45];
TR = TR(3);
t = TR.t;
Y = [-TR.LA; TR.TBy]';
Y = Y-min(Y);
Y = Y./max(abs(Y));
Y(:,2) = smooth(Y(:,2),91);
%Y(t<xlims(1) | t>xlims(2),:) = nan;

%event times
TT = [-0.1 0; 0 0.145; -0.070 0.225; .180 .305];

S = def_events({'LA clo';'LA rel';'PHAR [a]';'LA clo'},TT,[3 2 1 3]');

tvcolors = repmat([.5 .5 .5],2,1);
ecolors = N.color(ismember(N.label,{'clo' 'rel' '[a]'}),:);
S.color = ecolors([1 2 3 1],:);

%%
set(gcf,'currentaxes',ax(2));
h.he = draw_events(S,h);

k = 1e5;
beta = -2*sqrt(k);
tt = t;
    function dydt = td(t,y,F)
        f = F(find(tt>=t,1,'first'));
        dydt = [y(2); beta*y(2) + k*(f-y(1))];
    end

tspan = t;
y0 = [0; 0];
offr = 0;

for i=1:3
    yr = minmax(h.he{i}.fh.YData(:)');
  
    F = zeros(1,length(t));
    if i==1
        F(t>=TT(i,1) & t<=(TT(i,2)-offr)) = 1;
        F(t>=TT(end,1) & t<=(TT(end,2)-offr)) = 1;
    else
        F(t>=TT(i,1) & t<=(TT(i,2)-offr)) = 1;
    end
    
    [t,y] = ode45(@(t,y)td(t,y,F),tspan,y0);
    
    y = y*diff(yr);
    y = y+yr(1);
    plot(t,y(:,1),'linew',1,'color',ecolors(i,:));
    
    copyobj(h.he{i}.th,gca);
end

for i=1:2
    set(gcf,'currentaxes',ax(i+2));
    h.ph_tv(i) = plot(t,Y(:,i),'color',tvcolors(i,:),'linew',3); hold on;
    %h.th_tv(i) = text(xlims(1),Y(find(t>=xlims(1),1,'first'),i),tvlabs{i},'hori','right','fontsize',fs(2));
    h.th_tv(i) = text(xlims(1),Y(find(t>=xlims(1),1,'first'),i),tvlabs{i},'hori','left','verti','top','fontsize',fs(2));
    
end
axix = [3 3 4 3];

for i=1:size(TT,1)
    ix = t>=TT(i,1) & t<=TT(i,2);
    IX{i} = find(ix);
    T{i} = t(ix);
    h.lh_ev(i) = plot(T{i},Y(ix,axix(i)-2),'color',S.color(i,:),'linew',3,'parent',ax(axix(i)));
    h.fh_ev(i) = fill(TT(i,[1 2 2 1]),[0 0 1 1],S.color(i,:),'edgecolor','none','parent',ax(axix(i)),'facealpha',0.2);
end

axt = ax(2:end);
set(axt,'xlim',xlims);
axrescaley(0.05,axt);

set(gcf,'currentaxes',ax(2));
yo = 0.025; arlen = 0.3;
for i=1:size(TT,1)
    PP = [TT(i,1) max(h.he{i}.fh.YData)] + [0 arlen; 0 0.01];
    arrow(PP(1,:),PP(2,:),'length',8,'tipangle',30,'linewidth',3);
    %text(PP(1,1),PP(1,2)-yo,'?','fontsize',fs(1),'verti','bot','hori','center');
    PP = [TT(i,2) max(h.he{i}.fh.YData)] + [0 arlen; 0 0.01];
    arrow(PP(1,:),PP(2,:),'length',8,'tipangle',30,'linewidth',2);
    arrow(PP(1,:),PP(2,:),'length',8,'tipangle',30,'linewidth',0.5,'FaceColor','w');
    %text(PP(1,1),PP(1,2)-yo,'?','fontsize',fs(1),'verti','bot','hori','center');
end

%%
set(ax,'visible','off');
xlh = xlabel(ax(end),'time \rightarrow','fontsize',fs(1));
set(xlh,'visible','on');

stfig_panlab(ax(1:end-1),arrayfun(@(c){char(c+64)},1:length(ax)-1),'fontsize',fs(1));

%%

h.printfig(mfilename);


end

function [C] = draw_connections(C)
for i=1:height(C)
    draw_connection(C.conn{i}','linewidth',C.linewidth(i),'linestyle',C.linestyle{i},'conntype',C.conntype{i});
end

end

function [S] = def_events(labs,xpos,ypos)
S = table(labs,xpos,ypos,'VariableNames',{'label','xpos','ypos'});
end

%%
function [he] = draw_events(S,h)
hh = 0.65;
for i=1:height(S)
    he{i} = draw_event(S.xpos(i,:),S.color(i,:),S.label{i},'ypos',S.ypos(i),'height',hh,'endlines',false,'fontsize',h.fs(1));
end
end


%% base system
function [N] = base_system(h)

% node_info = {...
%     {'oris.','orb.'} [-1 0];
%     {'depr.','ang. oris'} [1 0];
%     'UL' [-1.5 1];
%     'LL' [0 1];
%     'JAW' [1 1];
%     'LA' [0 2];
%     'clo' [-1 3];
%     'rel' [1 3];
%     '/p/' [0 4]};

node_info = {...
    'UL' [0 1];
    'LL' [1 1];
    'JAW' [2 1];
    'TR' [3 1];
    'LA' [1 2];
    'PHAR' [2.5 2];
    'clo' [0.5 3];
    'rel' [1.5 3];
    '[a]' [2.5 3];
    '/p/' [1 4];
    '/a/' [2.5 4]};

N = table(node_info(:,1),cell2mat(node_info(:,2)),'VariableNames',{'label','pos'});

N.radius = 0.25*ones(height(N),1);
N.fontsize = h.fs(2)*ones(height(N),1);

%color by level
colors{4}  = lines(2);
colors{3}  = [pastelize(colors{4}(1,:),[0 .5]); colors{4}(2,:)];
colors{2}  = flipud(hsv(2));
colors{1} =  gray(4);

levels = unique(N.pos(:,2));
for i=1:length(levels)
    ixs = find(N.pos(:,2)==levels(i));
    N.color(ixs,:) = colors{i};
end

N.interpreter = repmat({'none'},height(N),1);

N.name = N.label;

for i=1:height(N)
    if iscell(N.name{i})
        N.name{i} = strjoin(N.name{i},'_');
    end
end

%connections:
% for i=1:height(N)
%     N.conn{i} = find(N.pos(:,2)==N.pos(i,2)-1);
% end
N = connect(N,'/p/',{'clo','rel'});
N = connect(N,'/a/',{'[a]'});
N = connect(N,'clo',{'LA'});
N = connect(N,'rel',{'LA'});
N = connect(N,'[a]',{'PHAR'});
N = connect(N,'LA',{'UL','LL','JAW'});
N = connect(N,'PHAR',{'JAW','TR'});

end

function [N] = connect(N,name1,name2)
N.conn{ismember(N.name,name1)} = find(ismember(N.name,name2));
end


function [h] = draw_event(xlims,varargin)

p = inputParser;

defaultcolor = 'r';
defaultlabel = ' ';
defaultypos = 0;
defaultheight = 1/3;
defaultfontsize = 36;
defaultedgecolor = 'none';
defaultfacealpha = 0.5;
defaultendlines = true;
defaultbegincolor = [.2 .6 .2];
defaultendcolor = 'r';
defaultendlinewidth = 4;
defaultendlineextent = 1.1;
defaultparent = gca;

addRequired(p,'xlims');
addOptional(p,'color',defaultcolor);
addOptional(p,'label',defaultlabel,@(x)ischar(x));
addParameter(p,'ypos',defaultypos);
addParameter(p,'height',defaultheight);
addParameter(p,'fontsize',defaultfontsize);
addParameter(p,'edgecolor',defaultedgecolor);
addParameter(p,'facealpha',defaultfacealpha);
addParameter(p,'endlines',defaultendlines);
addParameter(p,'begincolor',defaultbegincolor);
addParameter(p,'endcolor',defaultendcolor);
addParameter(p,'endlinewidth',defaultendlinewidth);
addParameter(p,'endlineextent',defaultendlineextent);
addParameter(p,'parent',defaultparent);


parse(p,xlims,varargin{:});
R = p.Results;

h.fh = fill(xlims([1 2 2 1]),R.ypos+[0 0 R.height*[1 1]],R.color,'FaceAlpha',R.facealpha,'EdgeColor',R.edgecolor); hold on;
%h.th = text(min(xlims)+0.01*mean(xlims),R.ypos+R.height/2,R.label,'verti','mid','fontsize',R.fontsize,'hori','left');
xo = 0.01*diff(get(R.parent,'xlim'));
h.th = text(min(xlims)+xo,R.ypos+R.height/2,R.label,'verti','mid','fontsize',R.fontsize,'hori','left');

h.allh = [h.fh h.th];

if R.endlines
    yf = R.height*(R.endlineextent-1)/2;
    ye = [-yf R.height+yf];
    h.elh(1) = line(xlims(1)*[1 1],R.ypos + ye,'color',R.begincolor,'linewidth',R.endlinewidth);
    h.elh(2) = line(xlims(2)*[1 1],R.ypos + ye,'color',R.endcolor,'linewidth',R.endlinewidth);
    h.allh = [h.allh h.elh];
end

end

