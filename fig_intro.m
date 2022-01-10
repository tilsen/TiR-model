function [] = fig_intro()

dbstop if error; close all;
h = fbmod_helpers;

ax = stf([1 1; 2 3],[0.03 -0.02 0.03 0.01],[0 0.05],'aspect',1);

%ax(2).Position(1) = ax(2).Position(1)-0.05;
ax(3).Position(1) = ax(3).Position(1)+0.05;

scales(ax(1),h);
units(ax(2),h);
systems(ax(3),h);

%%
stfig_panlab(ax,{'A','B','C'},'xoff',[0 0 0],'yoff',-[0.075 0.1 0.1],'hori','right');

%%
h.printfig(mfilename);

end

%% unit hierarchies and temporal projection
function [] = units(ax,h)

axes(ax);

U = {{'C' 'V' 'C'};...
    {'\mu' '\mu'}; ...
    {'\sigma' '\sigma'};...
    {'Ft' 'Ft'};...
    {'\omega' '\omega'};...
    {'iP' 'iP'};
    {'IP' 'IP'};
    {'Utt'}};

S = [];
c=1;
for i=1:size(U,1)
    ux = U{i};
    for j=1:length(ux)
        S(c).id = c;
        S(c).lab = ux{j};
        S(c).row = i;
        S(c).col = j;
        c=c+1;
    end
end

S = struct2table(S);

S.conn = nan(height(S),1);

S.conn(1:2) = find(ismember(S.lab,'\mu'),1,'first');
S.conn(3) = find(ismember(S.lab,'\mu'),1,'last');
for i=4:height(S)-1
    S.conn(i) = find(S.row==(S.row(i)+1),1,'last');
end

S.ypos = S.row-min(S.row);

xpos = arrayfun(@(c){linspace(-1,0,sum(S.row==c))},unique(S.row));

S.xpos = [xpos{:}]';

S.xpos = S.xpos.*(S.row.^1.5);

%correct CV
S.xpos(1:2) = S.xpos(4) + [-1 1];
S.ypos(1:3) = S.ypos(1:3) -0.5;

for i=6:height(S)
    if S.id(i)==S.id(find(S.row==S.row(i),1,'last'))
        S.xpos(i) = mean(S.xpos(S.row==(S.row(i)-1)));
    end
end

for i=1:height(S)-1
    xx = [S.xpos(i) S.xpos(S.conn(i))];
    yy = [S.ypos(i) S.ypos(S.conn(i))];
    plot(xx,yy,'k-','linew',2); hold on;
end

for i=1:height(S)
   S.th(i) = text(S.xpos(i),S.ypos(i),S.lab{i},'fontsize',h.fs(3),...
       'hori','center','verti','mid','BackGroundcolor','w',...
       'margin',0.01,'Fontname','cambria'); hold on; 
end

yo = -1;
pht = plot(xlim+[-1 1],yo*[1 1],'color','k','linew',2);
text(pht.XData(1),pht.YData(1),'time','verti','top','fontsize',h.fs(3));

ne = 5;
colors = flipud(gray(max(S.row)+ne));
colors = colors(ne:end,:);
S.color = colors(S.row,:);
S.lw = 0.5*S.row;

for i=1:height(S)
    S.ch{i} = children(S,S.id(i));
end

for i=6:height(S)
    
    if S.row(i)==1
        xx = sum(S.th(i).Extent([1 3]));
    elseif isempty(S.ch{i}) 
        xx = sum(S.th(i).Extent([1 3]));
    else
       xpos = S.xpos(S.ch{i}); 
       [~,ix] = max(xpos);
       Six = S.ch{i}(ix);
       xx = sum(S.th(Six).Extent([1 3]));
    end
    yy = S.th(i).Extent(2);
    plot(xx*[1 1],[yy yo-0.1*S.row(i)],'color',S.color(i,:),'linew',S.lw(i));
    
end

axis tight;
axrescaley(0.1);

set(gca,'Box','off','Visible','off');


end

function [ch] = children(S,i)

    ch = S.id(find(S.conn==i))';
    
    if ~isempty(ch)
        for j=1:length(ch)
            ch = [ch children(S,ch(j))];
        end
    end

end

%% systems schema
function [] = systems(ax,h)

axes(ax);

N.label = {'$x$'};
N.pos = [0 0];
N.radius = 1;
N.color = lines(1);

N(2).label = {''};
N(2).radius = 0.1;
N(2).pos = [1.5 1.5];
N(2).color = [1 1 1];

N = [N repmat(N(end),1,3)];
N = struct2table(N);
N.interpreter = repmat({'latex'},height(N),1);
N.fontsize = h.fs(1)*ones(height(N),1);

th = linspace(pi-pi/16,pi+pi/16,3)';
N.pos(3:end,:) = [1.5*cos(th) 1.5*sin(th)];

nh = draw_node(N);
axis equal;

C = [];
for i=2:length(nh)
    C = connect_objs(nh(i).fh,nh(1).fh);
    C = scale_connection(C',0.75)';
    draw_connection(C);
    delete(nh(i).fh);
end

text(-1.5,0,'$\mathcal{S}$','hori','right','fontsize',h.fs(2),'interpreter','latex');
text(1.5,1.5,'$Y$','hori','left','fontsize',h.fs(2),'interpreter','latex');

xp = 0;
yp = -2;
text(xp,yp,'change rule:','hori','center','fontsize',h.fs(2));
text(xp,yp-1.5,'$\frac{dx}{dt} = f(x) + F_\mathcal{S} + \sum_{j} F(y_j)$',...
    'hori','center','fontsize',h.fs(3),'interp','latex');


set(gca,'Visible','off','Clipping','off');
axis equal;
xlim(xlim+[-3 3]);
ylim(ylim+[-2 0]);

end

%% timescales
function [] = scales(ax,h)

axes(ax);

%%
theta_bounds = [3 8];
theta_per = 1./theta_bounds;
theta_range = theta_per/2;

fs = [36 18 16];

Nr = 11;

%%

tlims = [0.02 5];
tt = arrayfun(@(c){linspace(10^c,10^(c+1),10)},-2:0);
tt = unique([tt{:}]);
tt = tt(tt<=tlims(end));
plot(tlims,[0 0],'k-','linew',1); hold on;
oy = 0.25;
gcol = 0.75*ones(1,3);
for i=1:length(tt)
    plot(tt(i)*[1 1],oy*[-1 1],'k-');
    if tt(i)<0.1
        str = sprintf('%.02f',tt(i));
    elseif tt(i)<1
        str = sprintf('%0.1f',tt(i));
    else
        str = sprintf('%1.0f',tt(i));
    end
    str = strrep(str,'0.','.');
    th(i) = text(tt(i),-oy,str,...
        'verti','top','hori','center','fontsize',fs(3));
end

plot(theta_range(1).*[1; 1],oy*[-1; 1],'r-','linew',2);


set(th(ismember({th.String},{'.07','.70','7'})),'Visible','off');
set(th(ismember({th.String},{'.09','.90','9'})),'Visible','off');

set(gca,'XScale','log','Visible','off','ylim',[-Nr-1 Nr],'xlim',tlims);

th(2).String = {th(2).String, [blanks(20) 'timescale (sec)']};

%%

%---patterns
X(1,:) =     {'pat' 'movement velocity' [0.02 0.05]};
%X(end+1,:) = {'pat' 'movement range'    [0.02 0.10]};
X(end+1,:) = {'pat' 'CV relative timing'  [0.02 0.100]};
X(end+1,:) = {'pat' 'VC relative timing'  [0.100 0.350]};
X(end+1,:) = {'pat' 'segment durations'  [0.02 0.350]};
X(end+1,:) = {'pat' 'mora duration'     [0.10 0.250]};
X(end+1,:) = {'pat' 'syllable duration' [0.10 0.500]};
X(end+1,:) = {'pat' 'kinematic landmarks' [0.02 0.10]};
X(end+1,:) = {'pat' 'boundary-related duration' [0.5 5]};
X(end+1,:) = {'pat' 'speech rate (e.g. syllables/sec)' [1 4.5]};

%---vocabularies
X(end+1,:) = {'voc' 'intragestural coordination (stiffness)'  [0.02 0.10]};
X(end+1,:) = {'voc' 'intergestural coordination (CV)'  [0.02 0.100]};
X(end+1,:) = {'voc' 'intergestural coordination (VC)'  [0.100 0.350]};
X(end+1,:) = {'voc' 'gestural activation'  [0.02 0.250]};
X(end+1,:) = {'voc' 'phones'  [0.02 0.250]};
X(end+1,:) = {'voc' 'syllable constituents (onset, coda)'     [0.100 0.250]};
X(end+1,:) = {'voc' 'moras'     [0.100 0.250]};
X(end+1,:) = {'voc' 'syllables / syllable boundaries' [0.100 0.500]};
X(end+1,:) = {'voc' 'syllable boundaries' [0.100 0.500]};
%X(end+1,:) = {'voc' 'feet' [0.250 0.750]};
X(end+1,:) = {'voc' '\mu gestures (stress)' [0.300 0.750]};
X(end+1,:) = {'voc' 'coordination of \mu gestures' [0.300 0.750]};
X(end+1,:) = {'voc' 'prosodic words' [0.500 1.250]};
X(end+1,:) = {'voc' 'prosodic phrases / phrase boundaries' [1.0 3.0]};
X(end+1,:) = {'voc' '\pi gestures (phrase boundaries)' [1.0 3.0]};
X(end+1,:) = {'voc' 'coordination of \pi gestures' [1.0 3.0]};
X(end+1,:) = {'voc' 'intermediate phrases' [1.0 5.0]};
X(end+1,:) = {'voc' 'utterances' [1.0 5.0]};

%%
X = array2table(X,'VariableNames',{'type' 'label','t'});
X.t = cell2mat(X.t); 
X.t0 = X.t(:,1);
X.t1 = X.t(:,2);

types = unique(X.type,'stable');
offsets = [1 -Nr-1];
X.yo = cellfun(@(c)offsets(ismember(types,c)),X.type);

X = sortrows(X,{'type','t0', 't1'});

for i=1:length(types)
    ixs = ismember(X.type,types{i});
    X.type_ix(ixs) = (0:(sum(ixs)-1))';
end
X.y = X.yo + mod(X.type_ix,Nr-1);

X.color = repmat([gcol],height(X),1);
X.fs = fs(3)*ones(height(X),1);

hcolors = pastelize(flipud(lines(3)),0.25);

ix_boundaries = contains(X.label,'boundaries');
ix_coordination = contains(X.label,'coordination');
X.color(ix_boundaries,:) = repmat(hcolors(1,:),sum(ix_boundaries),1);
X.color(ix_coordination,:) = repmat(hcolors(2,:),sum(ix_coordination),1);

%%
tt = linspace(min(X.t(:)),max(xlim),10000);

for i=1:height(X)
    yp = X.y(i)+[0 1];
    xx = X.t(i,[1 2 2 1]);

    aa = ones(size(tt));
    tixs = [find(tt<=X.t(i,1),1,'last') find(tt>=X.t(i,2),1,'first')];
    
    sd0 = 0.1*X.t(i,1).^2;    
    tt0 = tt(1:tixs(1));
    aa0 = exp(-((tt0-tt(tixs(1))).^2)/sd0);
    aa0 = aa0/max(aa0);
    aa(1:tixs(1)) = aa0;
    
    sd1 = 0.1*X.t(i,2).^2;
    tt1 = tt(tixs(2):end);
    aa1 = exp(-((tt1-tt(tixs(2))).^2)/sd1);
    aa(tixs(2):end) = aa1;
    
    XX = [repmat(tt(1:end-1),2,1); repmat(tt(2:end),2,1)];
    YY = repmat(yp([1 2 2 1])',1,size(XX,2));
    
    patch(XX,YY,X.color(i,:),'EdgeColor','none','FaceVertexAlphaData',aa(1:end-1)','FaceAlpha','flat');
    X.xpos(i) = mean(xx);
    X.ypos(i) = X.y(i);
end

for i=1:height(X)
    th(i) = text(X.xpos(i),X.ypos(i),X.label{i},...
        'verti','bot','fontsize',X.fs(i),'hori','center');  

end

ix = find(contains({th.String},'rate'));
th(ix).Position(1) = th(ix).Position(1)-0.5;

ix = find(contains({th.String},'-related'));
th(ix).Position(1) = th(ix).Position(1)-0.5;

text(theta_range(1),oy,{'\phi=\pi,','{\it{f}}=3 Hz'},...
    'verti','bot','hori','center','fontsize',fs(3));

xoo = 0.001;
text(min(xlim)-xoo,max(ylim)-1.5,{'measurement','vocabulary'},...
    'verti','top','fontsize',fs(2),'fontangle','italic');
text(min(xlim)-xoo,mean(ylim)-3,{'theoretical','vocabulary'},...
    'verti','top','fontsize',fs(2),'fontangle','italic');


end


