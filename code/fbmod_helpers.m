function [h] = fbmod_helpers()

dbstop if error;

%add all project folders to path:
addpath(genpath('..'));

external_toolboxes = {
    'M:\Projects\toolboxes\Colormaps\'  'https://www.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps';
    'M:\Projects\toolboxes\arrow\'        'https://www.mathworks.com/matlabcentral/fileexchange/278-arrow), MATLAB Central File Exchange.' ;
    'M:\Projects\toolboxes\ds2nfu\'       'https://www.mathworks.com/matlabcentral/fileexchange/10656-data-space-to-figure-units-conversion';
    'M:\Projects\toolboxes\arrow\'        'https://www.mathworks.com/matlabcentral/fileexchange/278-arrow'};

%check for toolboxes
for i=1:size(external_toolboxes,1)
    if ~exist(external_toolboxes{i},'dir')
        fprintf('toolbox directory: %s not found. Model plotting functions will not work.\nObtain this toolbox from %s and update paths in this file\n',...
            external_toolboxes{i,1},external_toolboxes{i,2});
    else
        addpath(external_toolboxes{i});
    end
end

h.sim_dir = ['..' filesep 'data' filesep];
h.figures_dir = ['..' filesep 'figures' filesep];

%some default fontsizes
h.fs = [48 36 28 18];

%utility functions
h.latexify = @(str)cellfun(@(c){['${' c '}$']},str);
h.spacer = @(str)cellfun(@(c){['\ {' c '}\ ']},str);
h.fontsizer = @(str,fsz)cellfun(@(c){['{\fontsize{' num2str(fsz) '}' c '}']},str);
h.getix = @(T,name)find(ismember(T.name,name));
h.adjpos = @(M,name,adj)fbmod_adj_pos(M,name,adj);
h.printfig = @(c)printtif([h.figures_dir c]);

end

