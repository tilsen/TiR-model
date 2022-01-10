function [E] = format_cq_systems(E,varargin)

if iscell(E)
    for i=1:length(E)
        E{i} = format_cq_systems(E{i},varargin{:});
    end
    return;
end

p = inputParser;

def_radius = 0.1;
def_textloc = 'inside';
def_fontsize = 24;
def_colors = lines(height(E));
def_facealpha = 0.5;
def_interpreter = 'latex';

addRequired(p,'E');
addOptional(p,'radius',def_radius);
addOptional(p,'textloc',def_textloc);
addOptional(p,'fontsize',def_fontsize);
addOptional(p,'colors',def_colors);
addOptional(p,'facealpha',def_facealpha);
addOptional(p,'interpreter',def_interpreter);

parse(p,E,varargin{:});

Ne = height(E);

E.color = p.Results.colors;
E.radius = repmat(p.Results.radius,Ne,1);
E.textloc = repmat({p.Results.textloc},Ne,1);
E.fontsize = repmat(p.Results.fontsize,Ne,1);
E.facealpha = repmat(p.Results.facealpha,Ne,1);
E.interpreter = repmat({p.Results.interpreter},Ne,1);

E.label = E.names;

latexify = @(str)cellfun(@(c){['${' c '}$']},str);

switch(p.Results.interpreter)
    case 'latex'
        ix_latex =ismember(E.interpreter,'latex');
        if numel(ix_latex)>1
            E.label(ix_latex) = latexify(E.label(ix_latex));
        else
            E.label{ix_latex} = ['${' E.label{ix_latex}' '}$'];
        end
end


end