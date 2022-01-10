function [T] = fbmod_timer_constructor(timer_type,src,varargin)

p = inputParser;

def_trg = src;
def_thresh = nan;
def_valence = 1;
def_auto = false;
def_periodic = false;

addRequired(p,'timer_type');
addRequired(p,'src');
addParameter(p,'trg',def_trg);
addParameter(p,'thresh',def_thresh);
%addParameter(p,'valence',def_valence,@(x)isscalar(x) & ismember(x,[-1 0 1 nan]));
addParameter(p,'valence',def_valence,@(x)isscalar(x));
addParameter(p,'auto',def_auto);
addParameter(p,'periodic',def_periodic);

parse(p,timer_type,src,varargin{:});

T.type = {timer_type};
T.src = int16(src);
T.trg = int16(p.Results.trg);
T.thresh = p.Results.thresh;
T.valence = p.Results.valence;
T.auto = logical(p.Results.auto);
T.periodic = p.Results.periodic;

T = struct2table(T);

end