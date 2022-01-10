function [varargout] = stf(varargin)

dbstop if error;
p = inputParser;

if ~exist('inmarg','var'), inmarg=[]; end
if ~exist('exmarg','var'), exmarg=[]; end
if ~exist('axpan','var'), axpan=[]; end

default_axpan = [1 1];
default_aspect = Inf;
default_bare = true;
default_exmarg = [0.05 0.05 0.01 0.01];
default_inmarg = [0.05 0.05];
default_handlearray = 'vector';

addOptional(p,'axpan',default_axpan,@(x)isnumeric(x) || ishandle(x));
addOptional(p,'exmarg',default_exmarg,@(x)isnumeric(x));
addOptional(p,'inmarg',default_inmarg,@(x)isnumeric(x));
addParameter(p,'aspect',default_aspect);
addParameter(p,'bare',default_bare);
addParameter(p,'handlearray',default_handlearray);

%parse(p,axpan,exmarg,inmarg,varargin{:});
parse(p,varargin{:});

make_bare = @(h)set(h,'menubar','none','toolbar','none','name','','numbertitle','off');

if nargin==0
    figh = stfig_max();
    if p.Results.bare, make_bare(figh); end
    setaspect(p.Results.aspect,figh);
    varargout = {figh};
    return;
end

figh = stfig_max();
if p.Results.bare, make_bare(figh); end
setaspect(p.Results.aspect,figh);
ax = stfig_axpos(p.Results.axpan,[p.Results.exmarg p.Results.inmarg]);

switch(p.Results.handlearray)
    case 'matrix'
        if all(size(p.Results.axpan)==[1 2])
            ax = reshape(ax,p.Results.axpan(2),[])';
        else
            ax = ax(p.Results.axpan);
        end
end

varargout = {ax,figh};

end