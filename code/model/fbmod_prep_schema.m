function [M] = fbmod_prep_schema(M,varargin)

h = fbmod_helpers;

if iscell(M)
    for i=1:length(M)
        M{i} = fbmod_prep_schema(M{i},varargin{:});
    end
    return
end

p = inputParser;

%gesture, timer, sens
def_radius = [1/3 1/8 1/5];
def_fontsize = [24 20 20];
def_interpreter = {'latex' 'latex' 'latex'};
def_facealpha = [0.50 0.25 0.50];
def_labeltype = {'full' 'short' 'full'};
def_timertypes = {'itmr','etmr','stmr','extr','osc'};

def_offsets_Y       =  [0.0 -4]; %in G radii
def_offsets_T_extr  = [-1.0 4]; %in G radii
def_offsets_T_osc   =  [0.0 4]; %in G radii
def_offsets_T_etmr  =  [0.5 2]; %in G radii
def_offsets_T_stmr  =  [0.0 -2]; %in G radii

addRequired(p,'M');
addParameter(p,'radius',def_radius);
addParameter(p,'fontsize',def_fontsize);
addParameter(p,'interpreter',def_interpreter);
addParameter(p,'facealpha',def_facealpha);
addParameter(p,'labeltype',def_labeltype);
addParameter(p,'offsets_Y',def_offsets_Y);
addParameter(p,'offsets_T_extr',def_offsets_T_extr);
addParameter(p,'offsets_T_osc',def_offsets_T_osc);
addParameter(p,'offsets_T_etmr',def_offsets_T_etmr);
addParameter(p,'offsets_T_stmr',def_offsets_T_stmr);
addParameter(p,'timertypes',def_timertypes);

parse(p,M,varargin{:});

%expand in case of single element input
repf = {'radius' 'fontsize' 'interpreter' 'facealpha' 'labeltype'};
for i=1:length(repf)
    if numel(p.Results.(repf{i}))==1
        p.Results.(repf{i}) = repmat(p.Results.(repf{i}),1,3);
    end
end

% tables
M.G = make_gestures_table(M,p,h);
M.Y = make_sens_table(M,p,h);
M.T = make_timers_table(M,p,h);

%disable plotting of sensory systems which are unused:
M.Y.has_action = arrayfun(@(c)M.T.has_action(ismember(M.T.type,'stmr') & M.T.type_id==c),M.Y.type_id);


end

%% timers table
function [T] = make_timers_table(M,p,h)

T = M.S(ismember(M.S.type,p.Results.timertypes),:);
Nt = height(T);

T.src = T.id;
T.trg = arrayfun(@(c){find(M.chi(c,:)~=0)},T.src);

T.auto      = ismember(T.type,{'osc' 'extr'});
T.periodic  = ismember(T.type,'osc');
T.radius    = repmat(p.Results.radius(2),Nt,1);

T.G_ixs     = arrayfun(@(c)find(ismember(M.G.type_id,c)),T.type_id);
T.G_pos     = M.G.pos(T.G_ixs,:);
T.G_radius  = M.G.radius(T.G_ixs);
T.color     = M.G.color(T.G_ixs,:);

ix_itmr = ismember(T.type,'itmr');
ix_etmr = ismember(T.type,'etmr');
ix_stmr = ismember(T.type,'stmr');
ix_extr = ismember(T.type,'extr');
ix_osc = ismember(T.type,'osc');

T.pos(ix_itmr,:) = T.G_pos(ix_itmr,:) + T.G_radius(ix_itmr)*[0 1]  + T.radius(ix_itmr)*[0 -1];

T.pos(ix_etmr,:) = T.G_pos(ix_etmr,:) + p.Results.offsets_T_etmr.*T.G_radius(ix_etmr);
T.pos(ix_stmr,:) = T.G_pos(ix_stmr,:) + p.Results.offsets_T_stmr.*T.G_radius(ix_stmr);
T.pos(ix_osc,:)  = T.G_pos(ix_osc,:) + p.Results.offsets_T_osc.*T.G_radius(ix_osc);

%vertically stagger oscillators
T.pos(ix_osc,:)  = T.pos(ix_osc,:)  +  T.G_radius(ix_osc).*double([zeros(sum(ix_osc),1) mod(T.type_id(ix_osc),2)]);

%local extrinsic timers above all others
T.pos(ix_extr,:) = [T.G_pos(ix_extr,1) max(T.pos(:,2))*ones(sum(ix_extr),1)] + T.radius(ix_extr).*[-2 2];

switch(p.Results.labeltype{2})
    case 'short'
        T.label(ix_itmr) = arrayfun(@(c){['\widetilde{' num2str(c) '}']},T.type_id(ix_itmr));
        T.label(ix_etmr) = arrayfun(@(c){['\widehat{' num2str(c) '}']},T.type_id(ix_etmr));
        T.label(ix_stmr) = arrayfun(@(c){['\overline{' num2str(c) '}']},T.type_id(ix_stmr));
        T.label(ix_extr) = arrayfun(@(c){['\epsilon_{' num2str(c) '}']},T.type_id(ix_extr));
        T.label(ix_osc) = arrayfun(@(c){['\theta_{' num2str(c) '}']},T.type_id(ix_osc));
        
    case 'minimal'
        T.label(ix_itmr) = arrayfun(@(c){['~']},T.type_id(ix_itmr));
        T.label(ix_etmr) = arrayfun(@(c){['\wedge']},T.type_id(ix_etmr));
        T.label(ix_stmr) = arrayfun(@(c){['-']},T.type_id(ix_stmr));
        T.label(ix_extr) = arrayfun(@(c){['\epsilon']},T.type_id(ix_extr));
        T.label(ix_osc) = arrayfun(@(c){['\theta']},T.type_id(ix_osc));        
        
    otherwise
        T.label(ix_itmr) = arrayfun(@(c){['\widetilde{T}_{' num2str(c) '}']},T.type_id(ix_itmr));
        T.label(ix_etmr) = arrayfun(@(c){['\widehat{T}_{' num2str(c) '}']},T.type_id(ix_etmr));
        T.label(ix_stmr) = arrayfun(@(c){['\overline{T}_{' num2str(c) '}']},T.type_id(ix_stmr));
        T.label(ix_extr) = arrayfun(@(c){['\epsilon_{' num2str(c) '}']},T.type_id(ix_extr));
        T.label(ix_osc) = arrayfun(@(c){['\theta_{' num2str(c) '}']},T.type_id(ix_osc));
end

T.color(ix_extr,:) = repmat([.7 .7 .7],sum(ix_extr),1);

T.X = arrayfun(@(c){M.X(:,c)},T.id);

for i=1:height(T)
    if ~isempty(T.trg{i})
        T.F{i} = squeeze(M.F_XX(:,T.id(i),T.trg{i}));
        T.tau{i} = full(squeeze(M.tau(T.id(i),T.trg{i})));
        T.chi{i} = full(squeeze(M.chi(T.id(i),T.trg{i})));
    end
    switch(T.type{i})
        case 'osc'
            rix = M.getid(strrep(T.name{i},'osc','oscr1'));
            T.X{i} = M.X(:,rix).*cos(T.X{i});
    end
    
end

T.interpreter   = repmat(p.Results.interpreter(2),Nt,1);
T.fontsize      = repmat(p.Results.fontsize(2),Nt,1);
T.facealpha     = repmat(p.Results.facealpha(2),Nt,1);
T.style         = repmat({'node'},Nt,1);

T.style(ix_extr | ix_osc) = repmat({'box'},sum(ix_extr | ix_osc),1);
T.has_action = ~cellfun('isempty',T.F);

T.label(ismember(T.interpreter,'latex')) = h.latexify(T.label(ismember(T.interpreter,'latex')));

end

%% gestures table
function [G] = make_gestures_table(M,p,h)

G = M.S(ismember(M.S.type,'gest'),:);
Ng = height(G);

G.label = M.name';
G.color = M.color';

[ix_on,ix_off,t_on,t_off] = fbmod_gestural_activation(M.X(:,G.id),M.t);

G.ix_on = ix_on';
G.ix_off = ix_off';
G.t_on = t_on';
G.t_off = t_off';
G.t_on1 = cellfun(@(c)min([c nan]),G.t_on);
G.X = arrayfun(@(c){M.X(:,c)},G.id);

%default positions
G = sortrows(G,'t_on1');
G.pos = [(1:Ng)' zeros(Ng,1)];

G.interpreter   = repmat(p.Results.interpreter(1),Ng,1);
G.fontsize      = repmat(p.Results.fontsize(1),Ng,1);
G.radius        = repmat(p.Results.radius(1),Ng,1);
G.facealpha     = repmat(p.Results.facealpha(1),Ng,1);

G.label(ismember(G.interpreter,'latex')) = h.latexify(G.label(ismember(G.interpreter,'latex')));

end

%% sensory systems table
function [Y] = make_sens_table(M,p,~)

Y = M.S(ismember(M.S.type,'sens'),:);
Ny = height(Y);

Y.G_ix = arrayfun(@(c)find(ismember(M.G.type_id,c)),Y.type_id);

Y.label         = M.G.label(Y.G_ix);
Y.interpreter   = repmat(p.Results.interpreter(3),Ny,1);
Y.fontsize      = repmat(p.Results.fontsize(3),Ny,1);
Y.radius        = repmat(p.Results.radius(3),Ny,1);
Y.facealpha     = repmat(p.Results.facealpha(3),Ny,1);

Y.G_pos     = M.G.pos(Y.G_ix,:);
Y.G_radius  = M.G.radius(Y.G_ix);
Y.color     = M.G.color(Y.G_ix,:);

Y.pos = Y.G_pos + p.Results.offsets_Y.*Y.G_radius;

%Y.label(ismember(Y.interpreter,'latex')) = h.latexify(Y.label(ismember(Y.interpreter,'latex')));

end




% NG = sum(ismember(M.S.type,'gest'));
% G = make_gestures_table(NG);
%
%     function [G] = make_gestures_table(NG)
%        G = table({''},{''},nan,nan,nan,{''},{''},{[]},...
%            {[]},{[]},{[]},{[]},{''},nan,nan,nan,...
%            'VariableNames',{'name','type','type_id','id','x0','descr','label','color',...
%            'ix_on','ix_off','t_on','t_off','interpreter','fontsize','radius','facealpha'});
%        G = repmat(G,NG,1);
%     end
%
% G(:,[M.S.Properties.VariableNames]) = M.S(ismember(M.S.type,'gest'),:);
% G.label = M.name';
% G.color = M.color';
%
% [ix_on,ix_off,t_on,t_off] = fbmod_gestural_activation(M.X(:,G.id),M.t);
%
% G.ix_on = ix_on';
% G.ix_off = ix_off';
% G.t_on = t_on';
% G.t_off = t_off';
%
% NG = height(G);
%
% for i=1:NG
%     G.X{i} = X(:,G.id(i));
% end
%
% G.t_on1 = cellfun(@(c)min([c nan]),G.t_on);
%
% %default positions
% G = sortrows(G,'t_on1');
% G.pos = [(1:NG)' zeros(NG,1)];
%
% G.interpreter = repmat({p.Results.interpreter_G},NG,1);
% G.fontsize =    p.Results.fontsize_G*ones(NG,1);
% G.radius =      p.Results.radius_G*ones(NG,1);
% G.facealpha =   p.Results.facealpha_G*ones(NG,1);
%
% G.label(ismember(G.interpreter,'latex')) = h.latexify(G.label(ismember(G.interpreter,'latex')));

% %% gesture-specific sensory systems
% Y = G;
% Y.name = regexprep(Y.name,'gest_','sens_');
% Y.type = repmat({'sens'},height(Y),1);
% Y.id = cellfun(@(c)M.S.id(ismember(M.S.name,c)),Y.name);
% Y.type_id = cellfun(@(c)M.S.type_id(ismember(M.S.name,c)),Y.name);
%
% Y.interpreter = repmat({p.Results.interpreter_Y},height(Y),1);
% Y.fontsize = p.Results.fontsize_Y*ones(height(Y),1);
% Y.radius = p.Results.radius_Y*ones(height(Y),1);
% Y.facealpha = p.Results.facealpha_Y*ones(height(Y),1);
%
% g_ixs = arrayfun(@(c)find(ismember(G.type_id,c)),Y.type_id);
%
% Y.G_pos = G.pos(g_ixs,:);
% Y.G_radius = G.radius(g_ixs);
% Y.color = G.color(g_ixs,:);
%
% Y.pos = Y.pos - [zeros(height(Y),1) 4*Y.G_radius];