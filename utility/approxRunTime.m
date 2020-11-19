function t = approxRunTime(runtime, ncores, ntraj, popsize, niter, varargin)
% Estimate time (hrs) required for parameter optimisation, assuming time required
% depends mostly upon cost function evaluation, and given:
% runtime = time taken for single model run on single processor
% ncores = number of processors
% ntraj = number of trajectories (model runs)
% popsize = number of parameter sets used in population of optimiser
% niter = number of optimiser iterations

if ~isempty(varargin)
    idsp = strcmp(varargin, 'Display');
    if any(idsp), dsp = varargin{find(idsp)+1}; end
else
    dsp = true;
end

t = runtime * ntraj * (niter+1) * popsize / ncores / 60 / 60;

if dsp
    rt = round(t, 2, 'significant');
    disp(['optimsation run-time at least ' num2str(rt) ' hrs'])
end
