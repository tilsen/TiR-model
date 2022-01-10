function [C] = fbmod_connect_objs(H,M,varargin)

p = inputParser;

def_reduce_lines = [0.015 0.015]; %amount to reduce in normalized figure units
def_linewidth = 2;
def_connprops = {};

addRequired(p,'H');
addRequired(p,'M');
addOptional(p,'reduce_lines',def_reduce_lines);
addOptional(p,'linewidth',def_linewidth);
addOptional(p,'connprops',def_connprops);

parse(p,H,M,varargin{:});

rl = p.Results.reduce_lines;
connprops = [p.Results.connprops; {'linewidth', p.Results.linewidth}];

T = M.T;
[rdux,rduy] = nfu2axu(rl(1),rl(2));

C=[];
%connect etmr inputs from gestures
TX = T(ismember(T.type,{'etmr'}),:);
for i=1:height(TX)
    obj1 = H(ismember({H.name},TX.name{i})); 
    obj2 = H(ismember({H.name},sprintf('gest_%i',TX.type_id(i)))); 
    C = [C; connection_constructor(obj2,obj1,connprops{:})]; %#ok<*AGROW>
end

%connect stmr inputs from sensory systems
TX = T(ismember(T.type,{'stmr'}),:);
for i=1:height(TX)
    obj1 = H(ismember({H.name},TX.name{i})); 
    obj2 = H(ismember({H.name},sprintf('sens_%i',TX.type_id(i)))); 
    C = [C; connection_constructor(obj2,obj1,connprops{:})]; %#ok<*AGROW>
end

%connect oscillators if coupled
TX = T(ismember(T.type,{'osc'}) & T.has_action,:);
for i=1:height(TX)
    for j=1:height(TX)
        if M.phi(M.getid(TX.name{i}),M.getid(TX.name{j}))==0, continue; end
        obj1 = H(ismember({H.name},TX.name{i}));
        obj2 = H(ismember({H.name},TX.name{j}));
        C = [C; connection_constructor(obj1,obj2,connprops{:})]; %#ok<*AGROW>
    end
end
    
%generate table of all timer output connections
TX = T(ismember(T.type,{'etmr','stmr','extr','osc'}),:);
for i=1:height(TX)
    %lines to targets of actions
    trgs = TX.trg{i};
    obj1 = H(ismember({H.name},TX.name{i})); 
    for j=1:length(trgs)
        obj_name = M.S.name{ismember(M.S.id,trgs(j))};
        obj_name = strrep(obj_name,'gate_','gest_');
        obj_name = strrep(obj_name,'oscg_','osc_');
        obj_name = strrep(obj_name,'exg_','extr_');
        obj2 = H(ismember({H.name},obj_name)); 
        if isempty(obj2), continue; end
        C = [C; connection_constructor(obj1,obj2,connprops{:})];
    end
end

%generate lines and scale
for i=1:height(C)
    cc = connecting_line(C.obj1(i),C.obj2(i));
    cc_len = sqrt(sum(diff(cc).^2));
    rd_len = sqrt(sum(rdux^2+rduy^2));
    sc = (cc_len-rd_len)/cc_len;
    ccsc = scale_connection(cc,sc);
    C.conn{i} = ccsc;
end

end

%%
function [C] = connection_constructor(obj1,obj2,varargin)

if isempty(obj1) || isempty(obj2), C=[]; return; end    
    
p = inputParser;

def_linewidth = 2;
def_linecolor = [0 0 0];
def_conntype = 'forward';
def_label = '';

addRequired(p,'obj1');
addRequired(p,'obj2');
addParameter(p,'label',def_label);
addParameter(p,'linewidth',def_linewidth);
addParameter(p,'linecolor',def_linecolor);
addParameter(p,'conntype',def_conntype);

parse(p,obj1,obj2,varargin{:});

C.label = {p.Results.label};
C.obj1 = obj1.handles.fh;
C.obj2 = obj2.handles.fh;
C.name = {[obj1.name '_' obj2.name]};
C.name1 = {obj1.name};
C.name2 = {obj2.name};
C.linewidth = p.Results.linewidth;
C.linecolor = p.Results.linecolor;
C.conntype = {p.Results.conntype};


C = struct2table(C);

end



