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
% Set parFile = [] for default initial parameter values (hard-coded, not loaded)
% or parFile = 'filename.mat' to use saved values as the initials.
% parFile = [];
parFile = 'parameterInitialValues_RMS_Hellinger2_Atlantic_singleTraj_removeParams.mat';
Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', parFile);

display(Directories)

% Use default set-up to load and organise data and initialise model parameters.
% Set-up may be modified here by passing some extra arguments as name-value
% pairs (preferable), or directly modified within modelSetUp.m. % Useful 
% name-value pairs for modelSetUp.m include: useTraj, ESDmin, ESDmax, nsizes
% [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
%     'displayAllOutputs', true); % default set-up -- no plots
numTraj = 1; % use only a single trajectory per sampling event -- events at Arctic-Atlantic water boundaries are not easily descernable when using only one trajectory
[Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
    'displayAllOutputs', true, 'numTraj', numTraj);

% modelSetup.m may also produce plots -- set to 'true' any name-value pair as shown below
% [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
%     'plotCellConcSpectra', true, 'plotBioVolSpectra', true, ...
%     'plotSizeClassIntervals', true, ...
%     'trajectoryPlot', true, 'dendrogramPlot', true, ...
%     'plotScalarData', true, 'plotSizeData', true, 'plotAllData', true, ...
%     'displayForc', true, 'displayData', true);

% Parameters may be modified using name-value pairs in updateParameters.m
% Params = updateParameters(Params, FixedParams, 'pmax_a', 20, 'aP', 0.01);


%% Optimise parameters

% Choose which parameters to tune, the cost function, and numerical tuning
% algorithm. Can also select which data to fit / trajectories to run as
% Arctic and/or Atlantic.
niter = 50; % algorithm iterations (the process can be continued from where it stops, so this value doesn't matter that much)
popSize = 100; % number of parameter sets evaluated per iteration

costFunctionType = 'RMS_Hellinger2'; % fit scalar/nutrient data using least sum of errors, fit relative size data using Hellinger distances, totals in size data fit as ratio of zooplankton

% costFunctionType = 'RMSsmooth_Hellinger'; % fit scalar/nutrient data using least sum of squares (smoothed for robustness against model outliers), fit size data using Hellinger distances
% costFunctionType = 'meanCDFdist_Hellinger';
% costFunctionType = 'smoothCDFdist_Hellinger';
% costFunctionType = 'meanCDFdist_HellingerFullSpectrum';
% costFunctionType = 'meanCDFdist_HellingerFullSpectrum_averagedEventsDepths';

fitTrajectories = 'Atlantic';
% Different cost functions may use binned size data (integrated within
% modelled size class intervals) or may fit to the full size spectra data.
% This choice affects how model outputs are extracted to match data.
fitToFullSizeSpectra = false;
[FixedParams, Params, Forc, Data] = ...
    optimisationOptions(FixedParams, Params, Forc, Data, ...
    'niter', niter, ...
    'popSize', popSize, ...
    'costFunctionType', costFunctionType, ...
    'fitTrajectories', fitTrajectories, ...
    'fitToFullSizeSpectra', fitToFullSizeSpectra);

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

% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
v0 = initialiseVariables(FixedParams, Params, Forc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [size, depth, nutrient]
%                      organic matter      [type, depth, nutrient]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

restartRun = true; % restart algorithm from a saved prior run?
switch restartRun, case true
    fileName_results = 'fittedParameters';  % saved parameters file name
%     tag = '1';                              % and identifying tag
    tag = FixedParams.costFunction;
%     tag = [tag '_Atlantic_quadraticMortality_singleTraj_relativeSizeDataOnly'];
%     tag = [tag, '_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot_removeParams'];
    tag = [tag, '_Atlantic_singleTraj_removeParams'];
    fileName_results = fullfile(Directories.resultsDir, ...
        [fileName_results '_' tag]);
    % Load stored results    
    [~, results, ~, ~, boundsLower, boundsUpper, Data, Forc, FixedParams, Params, v0] = ...
        loadOptimisationRun(fileName_results);
    populationOld = results.populationHistory(:,:,end);
    scoresOld = results.scoreHistory(:,end);
    optimiserOptions = results.optimiserOptions;
    optimiserOptions.MaxGenerations = niter;
    optimiserOptions.InitialPopulationMatrix = populationOld;
    optimiserOptions.InitialScoresMatrix = scoresOld;
end

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Call optimiser
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = optimise(@(x) ... 
    costCalc(x, FixedParams, Params, Forc, Data, v0, FixedParams.odeIntegrator, ...
    FixedParams.odeOptions, 'selectFunction', costFunctionLabel), ... 
    npars, [], [], [], [], boundsLower, boundsUpper, [], optimiserOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
disp(['Optimisation time: ' num2str(floor(optimisationTime)) ' hrs, ' ...
    num2str(floor(mod(60*optimisationTime,60))) ' mins'])

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
    results.optPar = results_.optPar;
    results.optPar_summary = results_.optPar_summary;
    results.populationHistory = cat(3, results.populationHistory, ... 
        results_.populationHistory(:,:,2:end));
    results.scoreHistory = [results.scoreHistory, results_.scoreHistory(:,2:end)];
    results.optimiserOptions = results_.optimiserOptions;
else
    results = results_;
end

displayFittedParameters(results)

% Save output
saveParams = true;

fileName_results = 'fittedParameters';  % choose file name
tag = FixedParams.costFunction;         % and identifying tag
tag = [tag '_Atlantic_singleTraj_removeParams'];

tag = [tag '_Atlantic_singleTraj_adjustParBounds'];


fileName_results = fullfile(Directories.resultsDir, ...
    [fileName_results '_' tag]);

switch saveParams, case true
    % Fitted parameters are saved as the 'results' struct, and the 
    % associated model set-up is saved within each subsequent struct
    if ~exist('v0', 'var'), v0 = []; end
    saveOptimisationRun(fileName_results, results, Data, Forc, FixedParams, Params, v0);
end

% Fitted parameters may be saved as updated default values for future runs
updateStoredInitials = true;

fileName_initials = 'parameterInitialValues';
fileName_initials = fullfile(Directories.resultsDir, ... 
    [fileName_initials '_' tag]);

switch updateStoredInitials, case true
    saveFittedParams(fileName_initials, results.optPar, results.parNames)
end

