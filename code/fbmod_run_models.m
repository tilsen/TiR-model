function [M] = fbmod_run_models(M,varargin)

if ~iscell(M), M = {M}; end

p = inputParser;

def_parallel = false;
def_verbose = true;

addRequired(p,'M');
addOptional(p,'parallel',def_parallel);
addOptional(p,'verbose',def_verbose);

parse(p,M,varargin{:});

useParallel = p.Results.parallel;

if length(M)==1
    useParallel = false;
end

if useParallel
    poolobj = gcp('nocreate');
    if isempty(poolobj) 
        poolobj = gcp;
    end
    if poolobj.NumWorkers==1
        useParallel = false;
    end
end

L = length(M);

switch(useParallel)
    case true
        parfor i=1:L
            M{i} = fbmod_01(M{i});
        end
    case false
        for i=1:L
            status_str = status(sprintf('model %i/%i',i,L)); %#ok<NASGU>
            M{i} = fbmod_01(M{i});
        end
        status('clear');
end

if p.Results.verbose
    summarize_simulations(M);
end

if length(M)==1
    M = M{:};
end    

end

%%
function [] = summarize_simulations(M)

times = cellfun(@(c)c.elapsed_time,M);
iters = cellfun(@(c)c.Nt,M);
systems = cellfun(@(c)c.Ns,M);

hdrs = {'model simulations run:',...
    'avg. time (s) per simulation:',...
    'avg. time (ms) per iter:',...
    'avg. time (ms) per iter per system:'}';

hdrs = strjust(pad(hdrs),'right');

vals = [length(M) mean(times) 1000*mean(times./iters) 1000*mean(times./(iters.*systems))];

prec = {'%i','%1.1f','%1.1f','%1.3f'};

for i=1:length(hdrs)
    fprintf(['%s\t' prec{i} '\n'],hdrs{i},vals(i));
end

end