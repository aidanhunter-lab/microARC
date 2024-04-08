%% Size-structured 1D NPZD model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Numerically optimise model parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Set up model

% Refresh the workspace
clear; close all; delete(gcp('nocreate')); clc
% Include all subdirectories within search path
addpath(genpath(fileparts(which('fit_parameters'))))

rng(1) % set random seed

% Store folders/filenames of data and saved parameters.
% Set parFileName = [] for default initial parameter values (hard-coded, 
% not loaded) or parFileName = 'filename' to use saved values as the initials.
parFileName = [];
% parFileName = 'parameterInitialValues';
parFileType = '.mat';
% Identifying tag here is the name of cost function used to
% generate saved parameters, and a string describing model set-up
% (set tag=[] if unused).

tag = '_RMS_Hellinger_ZPratio_Atlantic_constraint';
% tag = '_RMS_Hellinger_ZPratio_Atlantic_correctedDiffusion';

% tag = '_RMS_Hellinger_ZPratio_Atlantic_final';
% tag = '_RMS_Hellinger2_Atlantic_aG_sigG_upweightAbnTot';

if ~isempty(parFileName), parFile = [parFileName tag parFileType]; else, parFile = []; end

Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', parFile);

display(Directories)

% Use default set-up to load and organise data and initialise model parameters.
% Set-up may be modified here by passing some extra arguments as name-value
% pairs (preferable), or directly modified within modelSetUp.m. % Useful 
% name-value pairs for modelSetUp.m include: useTraj, ESDmin, ESDmax, nsizes
numTraj = 1; % For optimisation efficiency use only a single trajectory per sampling event
[Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
    'displayAllOutputs', true, 'numTraj', numTraj);

% Parameters should be modified using name-value pairs in updateParameters.m
% and not by changing values directly in the Params struct, e.g.
% Params = updateParameters(Params, FixedParams, 'pmax_a', 20, 'aP', 0.01);


%% Optimise parameters
% Choose which parameters to tune, the cost function, and numerical tuning
% algorithm. Can also select which data to fit / trajectories to run as
% Arctic and/or Atlantic.

% Fit parameters to data sampled from either Atlantic or Arctic waters
fitTrajectories = 'Atlantic';
% List available cost function choices
[~,~,~,cfc] = costFunction();
% Select cost function from cfc
disp(cfc)
costFunctionType = 'RMS_Hellinger_ZPratio';
% Set rescaleForOptim true to estimate parameters within some transformed 
% space -- transforms chosen to improve optimisation efficiency by
% smoothing search space.
rescaleForOptim = true;
% Iterations of parameter-tuning algorithm (the algorithm can be halted
% then continued, so we can set niter quite small and run the algorithm in
% repeated blocks)
niter = 50;
% A population-based algorithm is used to tune parameters. Choose
% population size -- number of parameter sets evaluated per iteration
popSize = 100;

[FixedParams, Params, Forc, Data] = optimisationOptions(...
    FixedParams, Params, Forc, Data, ...
    'niter', niter, ...
    'popSize', popSize, ...
    'costFunctionType', costFunctionType, ...
    'fitTrajectories', fitTrajectories, ...
    'rescaleForOptim', rescaleForOptim);

% Optional arguments (e.g. 'niter') may be included as name-value pairs,
% otherwise default values are used.
% It is important to specify 'costFunctionType' as one of the options
% available in costFunction.m.
% The 'fitTrajectories' option is also important. It is used to select sets
% of trajectories to use within the parameter optimisation. Set 
% fitTrajectories = 'Atlantic' or 'Arctic' to use trajectories originating
% from either region and the associated fitting data, or set
% fitTrajectories = 'all' to fit to all data using using trajectories
% originating from Arctic or Atlantic. Omit fitTrajectories or set it empty
% ([]) to optimise using all trajectories and fitting to size data that is
% aggregated over all samples.

% Create functions defining parameter constraints (for use in the genetic
% algorithm, ga, function).
Params = parameterConstraints(FixedParams, Params);

% Are constraints satisfied by selected/default parameters?
pars = arrayfun(@(z) Params.(z{1}), FixedParams.tunePars);
constraint_vals = Params.constraints(pars);
constraints_satisfied = constraint_vals <= 0;
if(~all(constraints_satisfied))
    warning('Selected parameters do not satisfy constraints specified in Params.constraints')
end
disp(constraints_satisfied)

% How do previous fitted parameter sets compare to the constraints?
Params0 = Params; FixedParams0 = FixedParams;
prevFit = load('fittedParameters_RMS_Hellinger_ZPratio_Atlantic_correctedDiffusion.mat');
FixedParams.tunePars = prevFit.FixedParams.tunePars;
Params = updateParameters(Params, FixedParams, prevFit.fittingOutput.optPar);
Params = parameterConstraints(FixedParams, Params);

pars = arrayfun(@(z) Params.(z{1}), FixedParams.tunePars);
constraint_vals = Params.constraints(pars);
constraints_satisfied = constraint_vals <= 0;
if(~all(constraints_satisfied))
    warning('Selected parameters do not satisfy constraints specified in Params.constraints')
end
disp(constraints_satisfied)
% The optimal param set breaks the constraint, as we knew already. But what
% about the parameter-fitting history?
popfit = prevFit.fittingOutput.populationHistory;
popfit = permute(popfit, [2,1,3]);
popfit = popfit(:,:);
ntest = size(popfit, 2);
mtest = nan(ntest, 4);
for i = 1:ntest
    pars = popfit(:,i);
    Params = updateParameters(Params, FixedParams, prevFit.fittingOutput.optPar);
    constraint_vals = Params.constraints(pars);
    mtest(i,:) = constraint_vals;
end

cmet = mtest <= 0;

all(cmet, 1)
any(cmet(:,1))
any(cmet(:,2))
any(cmet(:,3))
any(cmet(:,4))
% Constraint 2 was always met (mu_max not too large), constraints 1 and 4
% were never met (mu_max optimum was always at too low volumes and mu_max
% was always too low at large sizes).

Params = Params0;
FixedParams = FixedParams0;




% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
v0 = initialiseVariables(FixedParams, Params, Forc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [size, depth, nutrient]
%                      organic matter      [type, depth, nutrient]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Restart algorithm from a saved prior run?
% restartRun = true;
restartRun = false;
switch restartRun, case true
    fileName_results = 'fittedParameters'; % saved parameters file name (and identifying tag should already be defined)
    fileType_results = '.mat';
    file_results = fullfile(Directories.resultsDir, ...
        [fileName_results tag fileType_results]);    
    fp = false; % fp = true => restart using full population from prior run
                % fp = false => restart with best param set from prior run and otherwise random population
    ni = true; % ni = true => use new v0 values generated (above) using loaded parameters
               % ni = false => use saved v0 values from previous optimisation run (v0 possibly generated using sub-optimal params)
    % Note: to avoid overwriting stucts Data, Forc, Params or
    % FixedParams generated above in model set-up, set to false the
    % optional name-value pair argument 'loadMainStructs' or individually
    % set as false 'loadData', 'loadForc', 'loadParams' or 'loadFixedParams'
    restartParamFit(file_results, optimiserOptions, 'fp', fp, 'ni', ni)
%     restartParamFit(file_results, optimiserOptions, 'fp', fp, 'ni', ni, ...
%         'loadData', false, 'loadForc', false, 'loadMainStructs', false)
end

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Call optimiser

% Non-linear constraints need to appear in argument after boundsUpper. This
% argument must be a function accepting X as an argument and returning
% vectors C and Ceq, where C(X) <= 0 and Ceq(X) = 0 define the constraints.
 % Store these constraints in the Params struct.

tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = optimise(@(x) ... 
    costCalc(x, FixedParams, Params, Forc, Data, v0, FixedParams.odeIntegrator, ...
    FixedParams.odeOptions, 'selectFunction', costFunctionType), ... 
    npars, [], [], [], [], boundsLower, boundsUpper, Params.constraints, [], ...
    optimiserOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
disp(['Optimisation time: ' num2str(floor(optimisationTime)) ' hrs, ' ...
    num2str(floor(mod(60*optimisationTime,60))) ' mins'])

% fun = @(x) costCalc(x, FixedParams, Params, Forc, Data, v0, FixedParams.odeIntegrator, ...
%     FixedParams.odeOptions, 'selectFunction', costFunctionType);
% [optPar, fval, exitflag, output, population, scores] = ga(fun, ... 
%     npars, [], [], [], [], boundsLower, boundsUpper, [], [], optimiserOptions);


% fun = @(x) costCalc(x, FixedParams, Params, Forc, Data, v0, FixedParams.odeIntegrator, ...
%     FixedParams.odeOptions, 'selectFunction', costFunctionType);
% nvars = npars;
% Aineq = [];
% bineq = [];
% Aeq = [];
% beq = [];
% lb = boundsLower;
% ub = boundsUpper;
% nonlcon = [];
% intcon = [];
% options = optimiserOptions;
% 
% ga(fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,options)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Examine and store results.

% If algorithm was terminated before reaching 'niter' iterations then the
% only output will be 'gapopulationhistory' and 'gacosthistory', which are 
% returned via gaStoreHistory.m
stoppedEarly = ~exist('fval', 'var');

% Store results
results_ = storeTunedOutput(FixedParams.optimiser, gapopulationhistory, ... 
    gacosthistory, optimiserOptions, FixedParams, 'stoppedEarly', stoppedEarly);

% Was optimisation continued from a prior run? If so, then append the fresh
% results into the loaded 'results' struct.
appendResults = exist('results', 'var');
if appendResults
    results.optPar_searchSpace = results_.optPar_searchSpace;
    results.optPar = results_.optPar;
    results.optPar_summary = results_.optPar_summary;
    results.populationHistory = cat(3, results.populationHistory, ... 
        results_.populationHistory(:,:,2:end));
    results.scoreHistory = [results.scoreHistory, results_.scoreHistory(:,2:end)];
    results.optimiserOptions = results_.optimiserOptions;
else
    results = results_;
end
clear results_

displayFittedParameters(results)

% Save output
saveParams = true;

% Choose file name
fileName_results = 'fittedParameters';
fileType_results = '.mat';
% and identifying tag
% tag = 'Atlantic_aG_sigG_upweightAbnTot';
tag = 'Atlantic_correctedDiffusion';
tag = ['_' FixedParams.costFunction '_' tag];

file_results = fullfile(Directories.resultsDir, ...
    [fileName_results tag fileType_results]);

switch saveParams, case true
    % Fitted parameters are saved as the 'results' struct, and the 
    % associated model set-up is saved within each subsequent struct
    if ~exist('v0', 'var'), v0 = []; end
    saveOptimisationRun(file_results, results, Data, Forc, FixedParams, Params, v0);
end

% Fitted parameters may be saved as updated default values for future runs
updateStoredInitials = true;
fileName_initials = 'parameterInitialValues';
fileType_initials = '.mat';
file_initials = fullfile(Directories.resultsDir, ... 
    [fileName_initials tag fileType_initials]);

switch updateStoredInitials, case true
    saveFittedParams(file_initials, results.optPar, results.parNames)
end

