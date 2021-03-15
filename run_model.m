%% Size-structured 1D NPZD model

clear
clc
close all
delete(gcp('nocreate'))

rng(1) % set random seed

%##########################################################################
%##########################################################################

%% Directories and files

% Choose forcing data
forcModel = 'sinmod';
% Forcing name (subfolder)
forcName = 'FramStrait_dep0-100';
% Experiment name
expName = 'AWI-Hausgarten';
% Forcing data directories
forcDir = fullfile('DATA', 'particle_trajectories', forcModel, forcName);
forcDummy = 'particles_MODEL_DOMAIN_EXPERIMENT-YEAR_YEAR_t*iSub03.mat';

% In-situ data directory
obsDir = fullfile('DATA', 'AWI_Hausgarten');

% Choose model type -- this will influence model set-up and some aspects of
% the ODEs.
bioModelOptions = {'singlePredatorClass', 'multiplePredatorClasses', 'mixotrophy'};
bioModel = bioModelOptions{1}; % Only singlePredatorClass is available so far

%##########################################################################
%##########################################################################

%% Model set-up

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Model input is stored in 4 structs: Forc, FixedParams, Params, and Data.
% Forc:        forcing data from physical model particle trajectories.
% FixedParams: constant model parameters (not numerically tuned).
% Params:      model parameters that MAY be tuned.
% Data:        observed data included in cost function to tune parameters.
% These structs are created in this code section.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select trajectories to use from physical model
% - options: 1. Use all trajectories by setting useTraj = []. Trajectories
%               are filtered later during model set-up.
%            2. Save RAM by manually defining trajectory index vector, 
%                e.g., useTraj = 1:50:5000 to extract every 50th trajectory
useTraj = [];
% useTraj = 1:200:5000;

% Specify years of available forcing data
years = 2017:2018;

% Load and extract relevant forcing data from physical model 'forcModel'.
F = loadForcing(forcModel, forcName, expName, years, forcDir, forcDummy, useTraj);

viewForcing = true;
switch viewForcing
    case true
        display(F); fprintf('\n\n')
        fields = fieldnames(F);
        disp('First year of forcing data. Dimensions include: yearday, depth layer, and particle trajectory.')
        display(F.(fields{1}))
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Load/prepare in-situ fitting data.
% prepareFittingData.m must be tailored to specific data sets.
Data = prepareFittingData(obsDir, ...
    'plotCellConcSpectra', true, 'plotBioVolSpectra', true);

% Choose number of modelled size class intervals using function
% chooseSizeClassIntervals.m
ESDmin = 1; % min/max sizes to retain from the data
ESDmax = 200;
nsizes = []; % number of modelled size classes (leave empty to view a range of values)
nsizesMin = 4; % min/max number of modelled size classes
nsizesMax = 12;
sizeData = chooseSizeClassIntervals(Data.size, 'datFull', Data.sizeFull, ...
    'ESDmin', ESDmin, 'ESDmax', ESDmax, ...
    'nsizes', nsizes, 'nsizesMin', nsizesMin, 'nsizesMax', nsizesMax, ...
    'plotSizeClassIntervals', true);
display(sizeData)

nsizes = 8; % number of modelled size classes (specifying a value makes chooseSizeClassIntervals.m return different output)
[Data.size, Data.sizeFull] = chooseSizeClassIntervals(Data.size, 'datFull', Data.sizeFull, ... 
    'ESDmin', ESDmin, 'ESDmax', ESDmax, ...
    'nsizes', nsizes, 'nsizesMin', nsizesMin, 'nsizesMax', nsizesMax, ...
    'plotSizeClassIntervals', true);

% Choose data types to use in cost function - I think bio-volume is the
% best choice for the size data
% Data = selectCostFunctionData(Data, ...
%     {'N','chl_a','PON','POC'}, {'BioVol'});
Data = selectCostFunctionData(Data, ...
    {'N','chl_a','PON','POC'}, {'BioVol'});

clear sizeData
close all

viewData = true; % display the fitting-data struct?
switch viewData
    case true
        display(Data); fprintf('\n\n')        
        disp('scalar data: Data.scalar')
        display(Data.scalar); fprintf('\n\n')
        disp('size-spectra (vector) data aggregated over sample events: Data.size')
        display(Data.size); fprintf('\n\n')
        disp('size-spectra data for separate sample events: Data.sizeFull')
        display(Data.sizeFull); fprintf('\n\n')
        disp('size-spectra data grouped into bins: Data.sizeFull.dataBinned')
        display(Data.sizeFull.dataBinned)
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Initialise model parameters.
% Values can be modified in defaultParameters.m, which is called from
% within initialiseParameters.m.
% Modelled cell size ranges are automatically chosen to correspond with the
% size class intervals already selected using the fitting data.
[FixedParams, Params] = initialiseParameters(F, Data);
display(FixedParams)
display(Params)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE:
% Parameters may be initialised (probably at better values) by loading
% saved values that resulted from optimisation, e.g.,

% [FixedParams, Params] = initialiseParameters(F, Data, ...
%     'load', 'results/parametersInitialValues_singlePredClass_2018');
% [FixedParams, Params] = initialiseParameters(F, Data, ...
%     'load', 'results/parametersInitialValues_singlePredClass_2018_2');


% [FixedParams, Params] = initialiseParameters(F, Data, ...
%     'load', ['results/parametersInitialValues_' bioModel]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Params0 = Params; % store the default values as separate workspace object

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE:
% Parameter values should not be changed by directly modifying Params
% because the size-dependent vector parameters will not update. Instead, 
% parameter values should be changed using name-value pair arguments in 
% updateParameters.m, e.g.,
% Params = updateParameters(Params, FixedParams, ... 
%     'pmax_a', 30, 'pmax_b', -0.55, 'k_G', 3);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Interpolate forcing data over modelled depth layers, and combine multiple
% years of forcing data into a single structure using prepareForcing.m.
% A few extra useful metrics are also calculated here.
Forc = prepareForcing(F,FixedParams);
clear F % to save memory

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% The 4 main structs (Forc, FixedParams, Params, Data) should now be stored
% in the workspace.
% The rest of this set-up code section selects which data to use, 
% standardises data, filters out unrequired data, and derives a few useful
% metrics.

% Remove fitting-data samples from below the maximum modelled depth
Data = omitDeepSamples(Data, FixedParams);

% Select which year(s) to model and filter out unused data.
% Function selectYears.m automatically chooses which year(s) to use based
% upon which years are most replete with data.
useSingleYear = true; % run model over a single year? If false, then 
                      % multiple years of forcing data MAY be used
                      % depending on fitting-data availability
useSingleSizeSpectra = true; % there are 2 size spectra for 2018, use only
                             % the Polarstern samples if true
[Data, Forc, FixedParams] = selectYears(Data, Forc, FixedParams, ... 
    'singleYear', useSingleYear, 'singleSpectra', useSingleSizeSpectra);

% There are lots of forcing data particle trajectories -- far too many to
% include within a parameter optimisation procedure, so these need to be
% filtered.
% For each in-situ sample event (specified in Data.scalar), choose a set of
% particle trajectories that pass nearby.
maxDist = 25; % Choose trajectories from within maxDist km radius of
              % sampling sites.
numTraj = 10; % Number of trajectories per sample event (duplicates may be 
              % required for some events, depending on maxDist).
[Forc_, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, numTraj);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE:
% The choice of maxDist is important. If maxDist is very large then the 
% trajectories ascribed to each sampling event will correspond poorly to
% the history of the water at the events, so small values are better.
% However, if maxDist is small then this can limit the number of
% trajectories available to use for each sampling event, and may entirely
% exclude some sampling events from the analyses if no trajectories are
% within the maxDist radius. We want to set maxDist small, but not too
% small...
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Observe any warnings produced by chooseTrajectories.m and use the table
% 'eventTraj' to help evaluate choice of maxDist and numTraj: if happy then
% set Forc = Forc_, otherwise reselect maxDist and numTraj and repeat
% chooseTrajectories.m
Forc = Forc_; clear Forc_

% Omit data from any sampling events not matched with a set of trajectories
[Data, eventTraj] = omitUnmatchedEvents(Data, eventTraj, Forc);

% For each trajectory, find the time of the latest sampling event.
% Integrations along trajectories can then be stopped at these events to
% reduce model run-times during parameter optimisation.
Forc = latestSampleTime(Forc, Data);

% Group sampling events by origin of particles: Arctic or Atlantic.
% Each individual trajectory is either of Atlantic or Arctic origin
% (stored in Forc.waterMass).
% Each sampling event is associated with a selection of trajectories, so
% may be considered as Atlantic, Arctic, or a mixture (stored in 
% Data.scalar.waterMass).
[Forc, Data.scalar] = particleOrigin(Forc, Data.scalar, ...
    'trajectoryPlot', true, 'dendrogram', true); pause(0.25)

% Standardise the fitting data using linear mixed models to adjust for
% variability due to depth and sampling event.
Data = standardiseFittingData(Data, ...
    'plotScalarData', true, 'plotSizeData', true, 'plotAllData', true);

close all

% Include extra fields indexing sorted order of data -- convenience for
% plotting
Data.scalar = sortOrderData(Data.scalar);

%##########################################################################
%##########################################################################

%% Integrate

% Initialise state variables.
% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
v0 = initialiseVariables(FixedParams, Params, Forc);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE:
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [depth, nutrient]
%                      organic matter      [type, depth, nutrient]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% % Use single precision -- get back to this... it's more involved...
% Forc = structDouble2Single(Forc);
% Params = structDouble2Single(Params);
% FixedParams = structDouble2Single(FixedParams);
% v0 = single(v0);

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj)
    poolObj = parpool('SpmdEnabled', false);
    numcores = poolObj.NumWorkers;
end

% Choose ODE solver.
integratorChoices = {'ode45', 'ode23', 'ode113', ...
    'ode15s', 'ode23s'};
odeIntegrator = integratorChoices{2};

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE:
% ode23 seems to be the most robust. ode45 occasionally produces NaNs.
% Maybe ode45 is less robust due to stiffness in the model equations. It
% seems that the equations are stiff for some parameter values, as the
% solver can slow down significantly... this is annoying, it would be
% useful to be able to dynamically switch between stiff and non-stiff
% solvers during the integrations (such switching functionality is
% available in the DEsolve package in R...)
% The stiff-solvers (ode15s, ode23s) are not working (at least partly)
% because they cannot use the non-negative constraint...
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Set solver options.
% Integrator is called separately for each time step -- defined by forcing 
% data (daily) intervals.
odeMaxTime = FixedParams.dt_max; % max timestep (days)
odeInitTime = 0.5 * odeMaxTime; % initial integration timestep (solver will
                                % automatically reduce this if required)

% odeOptions=odeset('AbsTol', 1e-6, 'RelTol', 1e-4,...
%     'InitialStep', odeInitTime, 'MaxStep', odeMaxTime, ...
%     'NonNegative', ones(FixedParams.nEquations,1));

odeOptions=odeset('AbsTol', 1e-6, 'RelTol', 1e-4,...
    'InitialStep', odeInitTime, 'MaxStep', odeMaxTime);

% Run the model
tic
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);
runtime = toc;
disp([num2str(runtime) ' seconds to integrate ' num2str(Forc.nTraj) ' trajectories']); fprintf('\n')
runtime = runtime / Forc.nTraj * numcores;
disp([num2str(runtime) ' seconds = average integration time of single trajectory on single processor']); fprintf('\n')

display(out)
display(auxVars)

%##########################################################################
%##########################################################################

%% Parameter tuning

% Select which parameters to optimise.
% Choose from the lists: Params.scalars & Params.sizeDependent.
parnames = {'wPOM', 'wp_a', 'wp_b', 'rDON', 'rPON', 'rPOC', 'beta2', ... 
    'beta3', 'aP', 'Gmax', 'k_G_a', 'k_G_b', 'pmax_a', 'pmax_b', ... 
    'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', 'Qmax_delQ_b', ... 
    'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b'};
npars = length(parnames);
par0 = nan(1,npars); % initial (default) values of tuning parameters
lb = nan(1,npars); % lower and upper bounds on tuning parameters
ub = nan(1,npars);
for i = 1:npars
    par0(i) = Params.(parnames{i});
    lb(i) = Params.bounds.(parnames{i})(1);
    ub(i) = Params.bounds.(parnames{i})(2);
end
FixedParams.tunePars = parnames;
FixedParams.tunePars_lb = lb;
FixedParams.tunePars_ub = ub;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select optimising algorithm (so far only ga is available)
optimiserChoices = {'ga','muga'};
optimiser = optimiserChoices{1};
optimise = str2func(optimiser);

% Select cost function.
% There's a few options for the cost function. Not yet sure which is the
% best... from this list choose either option 3 or 4
costFunctionChoices = { ...
    'LSS', ...
    'RMS', ...
    'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal', ...
    'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet', ...
    'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet' ...
    };
costFunctionType = costFunctionChoices{3}; % select cost function
FixedParams.costFunction = costFunctionType;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Test the cost function

% Parallelise integrations over trajectories
if isempty(gcp('nocreate'))
    poolObj = parpool('SpmdEnabled', false);
    numcores = poolObj.NumWorkers;
end

% Halt integrations at each trajectory's final sampling events to save time
Forc.integrateFullTrajectory = false;

% Integrate and calculate cost...
clear out auxVars
tic
[cost, costComponents, modData, out, auxVars] = costFunction(par0, ...
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, ... 
    'selectFunction', costFunctionType);
runtime = toc;
disp([num2str(runtime) ' seconds to integrate ' num2str(Forc.nTraj) ' trajectories']); fprintf('\n')
runtime = runtime / Forc.nTraj * numcores; 
disp([num2str(runtime) ' seconds = average integration time of single trajectory on single processor'])

display(modData) % modelled equivalents of fitting data
display(out) % state variable output
display(auxVars) % other model output
display(costComponents) % cost associated to each data type
display(cost) % total cost

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Optimise parameters using a numerical population-based algorithm

popSize = 100; % number of parameter sets in algorithm population
niter = 10; % algorithm iterations

% This approxRunTime.m function was supposed to estimate time required to
% run the optimisation, but it doesn't work because the ODE solver slows
% down when equations become stiff for some parameter values.
if exist('runtime', 'var'), approxOptimisationTime = ... 
        approxRunTime(runtime, numcores, Forc.nTraj, popSize, niter+1, 'Display', true); end

% Optimising algorithm options
gaOptions = optimoptions('ga', ...
    'PopulationSize', popSize, ...
    'Generations', niter, ...
    'InitialPopulationRange', [lb;ub], ...
    'InitialPopulationMatrix', par0, ...
    'Display', 'iter', ...
    'OutputFcn', @gaStoreHistory, ... % gaStoreHistory.m function dynamically stores output, which is available in workspace even if algorithm is halted prematurely
    'PlotFcn', @gaplotbestf);

% Call optimiser
clear cost costComponents modData out auxVars optPars fval
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = optimise(@(x) ... 
    costFunction(x, FixedParams, Params, Forc, Data, v0, odeIntegrator, ...
    odeOptions, 'selectFunction', costFunctionType), ... 
    npars, [], [], [], [], lb, ub, [], gaOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
disp(['Optimisation time: ' num2str(floor(optimisationTime)) ' hrs, ' ...
    num2str(floor(mod(60*optimisationTime,60))) ' mins'])

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Examine and store results.

% If algorithm was terminated before reaching 'niter' iterations then only
% 'gapopulationhistory' and 'gacosthistory' are returned via gaStoreHistory.m
stoppedEarly = ~exist('fval', 'var');

% Store results in gaOutput
switch stoppedEarly
    case false
        gaOutput = storeGAoutput(gapopulationhistory, gacosthistory, gaOptions, ...
            FixedParams, 'stoppedEarly', stoppedEarly);
    case true
        [gaOutput, optPar] = storeGAoutput(gapopulationhistory, gacosthistory, gaOptions, ...
            FixedParams, 'stoppedEarly', stoppedEarly);
end

viewResults = true;
switch viewResults
    case true
        display(gaOutput)
        if ~stoppedEarly
            bestFit = table(FixedParams.tunePars', FixedParams.tunePars_lb', ...
                optPar', FixedParams.tunePars_ub');
            bestFit.Properties.VariableNames = {'par','lower','opt','upper'};
            display(bestFit); fprintf('\n')
            disp(['Cost of best fit = ' num2str(fval)]); fprintf('\n')
            disp(['Optimisation exit flag = ' num2str(exitflag)]); fprintf('\n')
            display(output); fprintf('\n')
            display(population); fprintf('\n')
            display(scores)
        else
            display(gaOutput.optPar_summary); fprintf('\n')
            disp(['Cost of best fit = ' num2str(min(gaOutput.scoreHistory(:)))]); fprintf('\n')
        end
end

% Save output
tag = '_fitAllNutrientsAndBioVolPandZ';
fileName = ['results/fittedParameters_' FixedParams.costFunction tag];

saveParams = true;
switch saveParams
    case true
        saveOptimisationRun(fileName, gaOutput, Data, Forc, FixedParams, Params, v0);
end

% Should these parameter values be saved to use as initial values for
% future runs?
updateStoredInitials = false;
switch updateStoredInitials
    case true
        optPar = gaOutput.optPar;
        optPars_table = table(FixedParams.tunePars', optPar');
        optPars_table.Properties.VariableNames = {'Param','Value'};
        % take care choosing file name -- don't overwrite good values!
        saveObj = matfile(['results/parametersInitialValues_' bioModel]);
%         saveObj = matfile('results/parametersInitialValues_singlePredClass_2018_2');
        saveObj.Properties.Writable = true;
        saveObj.pars = optPars_table;
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% We may continue optimising parameters stored in gaOutput by initialising
% the algorithm with a stored population.

tag = '_fitAllNutrientsAndBioVolPandZ';
fileName = ['results/fittedParameters_' FixedParams.costFunction tag];

% Load stored results
[~, gaOutput, ~, ~, lb, ub, Data, Forc, FixedParams, Params, v0] = ... 
    loadOptimisationRun(fileName);


% Use stored result to initialise new optimisation run
populationOld = gaOutput.populationHistory(:,:,end);
scoresOld = gaOutput.scoreHistory(:,end);
gaOptions = gaOutput.gaOptions;

niter = 10;

gaOptions.MaxGenerations = niter;
gaOptions.InitialPopulationMatrix = populationOld;
gaOptions.InitialScoresMatrix = scoresOld;

clear optPars fval modData out auxVars
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = optimise(@(x) ...
    costFunction(x, FixedParams, Params, Forc, Data, v0, odeIntegrator, ... 
    odeOptions, 'selectFunction', costFunctionType), ... 
    npars, [], [], [], [], lb, ub, [], gaOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at');
disp(datetime('now')); fprintf('\n')
disp(['Optimisation time: ' num2str(floor(optimisationTime)) ' hrs, ' ...
    num2str(floor(mod(60*optimisationTime,60))) ' mins'])

% As above, examine and save output...
stoppedEarly = ~exist('fval', 'var');
switch stoppedEarly
    case false
        gaOutput = storeGAoutput(gapopulationhistory, gacosthistory, gaOptions, ...
            FixedParams, 'stoppedEarly', stoppedEarly);
    case true
        [gaOutput, optPar] = storeGAoutput(gapopulationhistory, gacosthistory, gaOptions, ...
            FixedParams, 'stoppedEarly', stoppedEarly);
end

viewResults = false;
switch viewResults
    case true
        display(gaOutput)
        if ~stoppedEarly
            bestFit = table(FixedParams.tunePars', FixedParams.tunePars_lb', ...
                optPar', FixedParams.tunePars_ub');
            bestFit.Properties.VariableNames = {'par','lower','opt','upper'};
            display(bestFit); fprintf('\n')
            disp(['Cost of best fit = ' num2str(fval)]); fprintf('\n')
            disp(['Optimisation exit flag = ' num2str(exitflag)]); fprintf('\n')
            display(output); fprintf('\n')
            display(population); fprintf('\n')
            display(scores)
        else
            display(gaOutput.optPar_summary); fprintf('\n')
            disp(['Cost of best fit = ' num2str(min(gaOutput.scoreHistory(:)))]); fprintf('\n')
        end
end

% Save output
saveParams = true;
switch saveParams
    case true
        saveOptimisationRun(fileName, gaOutput, Data, Forc, FixedParams, Params, v0);
end

updateStoredInitials = false;
switch updateStoredInitials
    case true
        optPar = gaOutput.optPar;
        optPars_table = table(FixedParams.tunePars', optPar');
        optPars_table.Properties.VariableNames = {'Param','Value'};
        % take care choosing file name -- don't overwrite good values!        
        saveObj = matfile(['results/parametersInitialValues_' bioModel]);
%         saveObj = matfile('results/parametersInitialValues_singlePredClass_2018_2');
        saveObj.Properties.Writable = true;
        saveObj.pars = optPars_table;
end

%##########################################################################
%##########################################################################

%% Plots

% Code folding for switches is not enabled by default, so type 'preferences'
% into command window then go to MATLAB -> Editor/Debugger -> Code Folding
% to enable folding on switches. This is useful but not necessary...
preferences('Editor/Debugger')

% Directory for saved plots
folder = 'results/plots/';

% Load fitted parameters, and associated data and fixed parameters...
tag = '_fitAllNutrientsAndBioVolPandZ';
fileName = ['results/fittedParameters_' FixedParams.costFunction tag];

[~, gaOutput, parnames, optPar, lb, ub, Data, Forc, FixedParams, Params, v0] = ... 
    loadOptimisationRun(fileName);


% Params = updateParameters(Params, FixedParams, ...
%     'Gmax', 5, 'k_G_a', 0.5, 'k_G_b', 0.18);
% optPar = cellfun(@(x) Params.(x), FixedParams.tunePars);

% Generate model output using fitted parameters
Params = updateParameters(Params, FixedParams, optPar);
costFunctionType = FixedParams.costFunction;
Forc.integrateFullTrajectory = true;

[cost, costComponents, modData, out, auxVars] = costFunction(optPar, ...
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, ... 
    'selectFunction', costFunctionType, 'returnExtra', 'all');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE:
% All of these plotting functions require further work to make them more
% generally applicable for any model outputs...
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~~~~~~~~~
% Model fit to data
%~~~~~~~~~~~~~~~~~~

save = true;

% Summary plots displaying model fit to data
logPlot = true; % for scalar data choose logPlot = true or false
pltChl = plot_fitToData('chl_a', Data, modData, logPlot); pause(0.25)
pltPON = plot_fitToData('PON', Data, modData, logPlot); pause(0.25)
pltPOC = plot_fitToData('POC', Data, modData, logPlot); pause(0.25)
logPlot = false;
pltN = plot_fitToData('N', Data, modData, logPlot); pause(0.25)

logPlot = 'loglog'; % for size spectra data choose logPlot = 'loglog' or 'semilogx'
pltCellConc = plot_fitToData('CellConc', Data, modData, logPlot); pause(0.25)
pltBioVol = plot_fitToData('BioVol', Data, modData, logPlot); pause(0.25)
pltNConc = plot_fitToData('NConc', Data, modData, logPlot); pause(0.25)

switch save, case true
    % scalar data
    if exist('pltChl', 'var') && isvalid(pltChl)
        filename = 'fitToData_chl.png';
        print(pltChl, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltN', 'var') && isvalid(pltN)
        filename = 'fitToData_DIN.png';
        print(pltN, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltPON', 'var') && isvalid(pltPON)
        filename = 'fitToData_PON.png';
        print(pltPON, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltPOC', 'var') && isvalid(pltPOC)
        filename = 'fitToData_POC.png';
        print(pltPOC, fullfile(folder, filename), '-r300', '-dpng');
    end
    
    % size data
    if exist('pltCellConc', 'var') && isvalid(pltCellConc)
        filename = 'fitToData_CellConcSpectra.png';
        print(pltCellConc, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol', 'var') && isvalid(pltBioVol)
        filename = 'fitToData_BioVolSpectra.png';
        print(pltBioVol, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltNConc', 'var') && isvalid(pltNConc)
        filename = 'fitToData_NConcSpectra.png';
        print(pltNConc, fullfile(folder, filename), '-r300', '-dpng');
    end
end


%~~~~~~~~~~~~~~~~~~
% Fitted parameters
%~~~~~~~~~~~~~~~~~~

save = false;

% Display fitted parameters in relation to their bounding values (in the
% table, columns widths shoukd be adjustable).
plt = plot_fittedParameters(gaOutput.optPar_summary);

switch save, case true
    if exist('plt', 'var') && isvalid(plt)
        filename = 'fittedParameters.png';
        print(plt, fullfile(folder, filename), '-r300', '-dpng');
    end
end


close all

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Contour plots -- single trajectories, or grouped by sample event
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

save = false;

% Choose one or more trajectories -- if multiple are selected then the plot
% will average over them.
sampleEvent = 10;
% All trajectories used for sampleEvent
traj = find(Data.scalar.EventTraj(sampleEvent,:));
% If waterMass is either Atlantic OR Arctic then it may make sense to plot
% average over all trajectories, although there could be unwanted smoothing
% effects...
% Otherwise, if waterMass is a mixture of origins, then group trajectories
% by origin and make separate plots
waterMass = Data.scalar.waterMass{sampleEvent};

switch waterMass
    case {'Atlantic', 'Arctic'}
        plt_Forc = plot_contour_DepthTime('forcing', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_DIN = plot_contour_DepthTime('DIN', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_OM = plot_contour_DepthTime('DOM_POM', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_N = plot_contour_DepthTime('phytoplankton_N', ...
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_Chl = plot_contour_DepthTime('phytoplankton_Chl', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_C = plot_contour_DepthTime('phytoplankton_C', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_N_C = plot_contour_DepthTime('phytoplankton_N_C', ...
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_Chl_N = plot_contour_DepthTime('phytoplankton_Chl_N', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_Z = plot_contour_DepthTime('zooplankton', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
    case 'Arctic/Atlantic'
        if strcmp(waterMass, 'Arctic/Atlantic')
            trajAtlantic = traj(strcmp(Forc.waterMass(traj), 'Atlantic'));
            trajArctic = traj(strcmp(Forc.waterMass(traj), 'Arctic'));
        end
        plt_Forc_Atl = plot_contour_DepthTime('forcing', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ... 
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_Forc_Arc = plot_contour_DepthTime('forcing', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ... 
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_DIN_Atl = plot_contour_DepthTime('DIN', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_DIN_Arc = plot_contour_DepthTime('DIN', ... 
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_OM_Atl = plot_contour_DepthTime('DOM_POM', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_OM_Arc = plot_contour_DepthTime('DOM_POM', ... 
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_N_Atl = plot_contour_DepthTime('phytoplankton_N', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_N_Arc = plot_contour_DepthTime('phytoplankton_N', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_Chl_Atl = plot_contour_DepthTime('phytoplankton_Chl', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_Chl_Arc = plot_contour_DepthTime('phytoplankton_Chl', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_C_Atl = plot_contour_DepthTime('phytoplankton_C', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_C_Arc = plot_contour_DepthTime('phytoplankton_C', ... 
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_N_C_Atl = plot_contour_DepthTime('phytoplankton_N_C', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_N_C_Arc = plot_contour_DepthTime('phytoplankton_N_C', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_Chl_N_Atl = plot_contour_DepthTime('phytoplankton_Chl_N', ...
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_Chl_N_Arc = plot_contour_DepthTime('phytoplankton_Chl_N', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_Z_Atl = plot_contour_DepthTime('zooplankton', ...
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_Z_Arc = plot_contour_DepthTime('zooplankton', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
end

switch save, case true
    
    switch waterMass
        
        case {'Atlantic', 'Arctic'}
        
            if exist('plt_Forc', 'var') && isvalid(plt_Forc)
                filename = ['forcing_data_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_Forc, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_DIN', 'var') && isvalid(plt_DIN)
                filename = ['DIN_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_DIN, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_OM', 'var') && isvalid(plt_OM)
                filename = ['OM_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_OM, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_N', 'var') && isvalid(plt_P_N)
                filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_N, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_Chl', 'var') && isvalid(plt_P_Chl)
                filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_Chl, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_C', 'var') && isvalid(plt_P_C)
                filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_C, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_N_C', 'var') && isvalid(plt_P_N_C)
                filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_N_C, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_Chl_N', 'var') && isvalid(plt_P_Chl)
                filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_Chl_N, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_Z', 'var') && isvalid(plt_Z)
                filename = ['zooplankton_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_Z, fullfile(folder, filename), '-r300', '-dpng');
            end
            
        case  'Arctic/Atlantic'
            
            if exist('plt_Forc_Atl', 'var') && isvalid(plt_Forc_Atl)
                filename = ['forcing_data_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_Forc_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Forc_Arc', 'var') && isvalid(plt_Forc_Arc)
                filename = ['forcing_data_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_Forc_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end

            if exist('plt_DIN_Atl', 'var') && isvalid(plt_DIN_Atl)
                filename = ['DIN_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_DIN_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_DIN_Arc', 'var') && isvalid(plt_DIN_Arc)
                filename = ['DIN_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_DIN_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end

            if exist('plt_OM_Atl', 'var') && isvalid(plt_OM_Atl)
                filename = ['OM_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_OM_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_OM_Arc', 'var') && isvalid(plt_OM_Arc)
                filename = ['OM_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_OM_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_N_Atl', 'var') && isvalid(plt_P_N_Atl)
                filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_N_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_N_Arc', 'var') && isvalid(plt_P_N_Arc)
                filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_N_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_Chl_Atl', 'var') && isvalid(plt_P_Chl_Atl)
                filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_Chl_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_Chl_Arc', 'var') && isvalid(plt_P_Chl_Arc)
                filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_Chl_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_C_Atl', 'var') && isvalid(plt_P_C_Atl)
                filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_C_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_C_Arc', 'var') && isvalid(plt_P_C_Arc)
                filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_C_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_N_C_Atl', 'var') && isvalid(plt_P_N_C_Atl)
                filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_N_C_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_N_C_Arc', 'var') && isvalid(plt_P_N_C_Arc)
                filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_N_C_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_Chl_N_Atl', 'var') && isvalid(plt_P_Chl_N_Atl)
                filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_Chl_N_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_Chl_N_Arc', 'var') && isvalid(plt_P_Chl_N_Arc)
                filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_Chl_N_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_Z_Atl', 'var') && isvalid(plt_Z_Atl)
                filename = ['zooplankton_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_Z_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Z_Arc', 'var') && isvalid(plt_Z_Arc)
                filename = ['zooplankton_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_Z_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
    end
end
    

close all

%~~~~~~~~~~~~~~~
% Time evolution
%~~~~~~~~~~~~~~~

% These polygon plots should be extended to show water masses of Atlantic
% and of Arctic orgin because some sampling events ue trajectories
% originating from both oceans... red for Atlantic blue for Arctic

save = false;

% Choose event
sampleEvent = 22;
if ~ismember(sampleEvent, 1:Data.scalar.nEvents), warning(['Choose event number within range (1, ' num2str(Data.scalar.nEvents) ')']); end
% trajectory indices
traj = find(Data.scalar.EventTraj(sampleEvent,:));

highlightColour = [1 0 1];
plotOptions = {'forcing', 'DIN', 'organicN', 'organicC', 'phytoplankton_C', ...
    'phytoplanktonStacked', 'phytoZooPlanktonStacked'};

for varIndex = 1:length(plotOptions)
    
    plotVariables = plotOptions{varIndex};
    
    switch plotVariables
        
        case 'forcing'
            % Forcing data
            plt_forcing = figure;
            plt_forcing.Units = 'inches';
            plt_forcing.Position = [0 0 12 12];
            
            % temperature
            subplot(3,1,1)
            plot_timeSeries_trajectoryPolygon('temperature', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            % diffusivity
            subplot(3,1,2)
            plot_timeSeries_trajectoryPolygon('diffusivity', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            % PAR
            subplot(3,1,3)
            plot_timeSeries_trajectoryPolygon('PAR', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'DIN'
            % DIN
            plt_DIN = figure;            
            plt_DIN.Units = 'inches';
            plt_DIN.Position = [0 0 12 8];
            
            subplot(2,1,1)
            plot_timeSeries_trajectoryPolygon('DIN', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,1,2)
            plot_timeSeries_trajectoryPolygon('DIN', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'organicN'
            plt_ON = figure;
            plt_ON.Units = 'inches';
            plt_ON.Position = [0 0 24 8];
            
            subplot(2,2,1)
            plot_timeSeries_trajectoryPolygon('DON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,2)
            plot_timeSeries_trajectoryPolygon('PON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,3)
            plot_timeSeries_trajectoryPolygon('DON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,4)
            plot_timeSeries_trajectoryPolygon('PON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'organicC'
            plt_OC = figure;
            plt_OC.Units = 'inches';
            plt_OC.Position = [0 0 24 8];

            subplot(2,2,1)
            plot_timeSeries_trajectoryPolygon('DOC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,2)
            plot_timeSeries_trajectoryPolygon('POC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,3)
            plot_timeSeries_trajectoryPolygon('DOC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,4)
            plot_timeSeries_trajectoryPolygon('POC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'phytoplankton_C'
            plt_P_C = plot_timeSeries_trajectoryPolygon('phytoplankton_C', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'max', 'highlightColour', highlightColour, 'fixedYaxis', false);
            plt_P_C_fixedScale = plot_timeSeries_trajectoryPolygon('phytoplankton_C', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'max', 'highlightColour', highlightColour, 'fixedYaxis', true);
            
        case 'phytoplanktonStacked'
            plt_P_C_stacked = plot_timeSeries_trajectoryPolygon('phytoplanktonStacked', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data);

        case 'phytoZooPlanktonStacked'
            plt_P_C_stacked = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data);

    end
    
end

switch save
    
    case true
        
        if exist('plt_forcing', 'var') && isvalid(plt_forcing)
            filename = ['forcing_data_timeSeriesPolygon_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_forcing, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_DIN', 'var') && isvalid(plt_DIN)
            filename = ['DIN_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_DIN, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_OM', 'var') && isvalid(plt_OM)
            filename = ['OM_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_OM, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_N', 'var') && isvalid(plt_P_N)
            filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_N, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_Chl', 'var') && isvalid(plt_P_Chl)
            filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_Chl, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_C', 'var') && isvalid(plt_P_C)
            filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_C, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_N_C', 'var') && isvalid(plt_P_N_C)
            filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_N_C, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_Chl_N', 'var') && isvalid(plt_P_Chl_N)
            filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_Chl_N, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_Z', 'var') && isvalid(plt_Z)
            filename = ['zooplankton_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_Z, fullfile(folder, filename), '-r300', '-dpng');
        end
end


plt = plot_timeSeries_barplot('phytoZooPlankton',sampleEvent,traj,out,FixedParams,Forc,Data);

switch save

    case true
        
        if exist('plt', 'var') && isvalid(plt)
            filename = ['barplot_timeSeries_plankton_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt, fullfile(folder, filename), '-r300', '-dpng');
        end
end


close all

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot groups of trajectories corresponding to each event
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% These plots need improved... better to run model over ALL available
% trajectories (not just those used for optimisation), then display output
% given that greater spatial coverage, and optionally show sampling events
% Also, move the plot functions to their own separate script to get rid of
% the outputPlot.m that aggregates them all...

% Choose event
ie = 7;
if ~ismember(ie, 1:Data.scalar.nEvents), warning(['Choose event number within range (1, ' num2str(Data.scalar.nEvents) ')']); end
% trajectory indices
kk = find(Data.scalar.EventTraj(ie,:));

outputPlot('trajectoryLine_LatLong','direction',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','forcing',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','inorganicNutrient',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','DOM_POM',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','phytoplankton_N',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','zooplankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);


close all

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% %### in progress...
% 
% ie = 1; % sampling event
% kk = find(Data.scalar.EventTraj(ie,:)); % all trajectories for selected event
% k = kk(7); % single trajectory
% 
% % x = out.N(:,:,:,:);
% x = out.N(:,:,:,k);
% 
% xmin = squeeze(min(x));
% xmax = squeeze(max(x));
% % tt = yearday(Forc.t(:,:));
% tt = yearday(Forc.t(:,k));
% n = size(tt, 1);
% 
% lat = Forc.y(:,k); % latitude
% % lat = Forc.y(:,:); % latitude
% % plot colour changing with latitude along line..?
% 
% xmin = reshape(xmin, [1 n]);
% xmax = reshape(xmax, [1 n]);
% tt = reshape(tt, [1 n]);
% lat = reshape(lat, [1 n]);
% z = zeros(size(xmin));
% 
% figure
% surface([tt;tt], [xmin;xmin], [z;z], [lat;lat], ...
%     'facecol', 'no', 'edgecol', 'interp', 'linew', 2)
% hold on
% surface([tt;tt], [xmax;xmax], [z;z], [lat;lat], ...
%     'facecol', 'no', 'edgecol', 'interp', 'linew', 2)
% hold off
% cb = colorbar();
% deg = char(176);
% cb.Label.String = ['latitude (' deg 'N)'];
% xlabel('day of year')
% ylabel(['DIN (mmol N m^{-3})'])
% title('max and min DIN over time')
% 
% 
% colormap(plasma)
% 
% % ###
% % Aggregate trajectories originating from Atlantic and Arctic waters then
% % make summmary plots displaying typical (across-trajectories) results...
% 
% pltN = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, 'DIN', 'minmax');
% 
% pltOM_N = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, 'OM_N', 'max');
% pltOM_C = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, 'OM_C', 'max');
% 
% now do the zooplanton... save phyto till last...
% 

