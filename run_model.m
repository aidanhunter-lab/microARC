%% Size-structured 1D NPZD model with a single predator class

clear
clc
close all

% rng(1) % set random seed

%% Load and filter forcing data

% select your model names
% bioModel  = 'size_based_model'; % What biogeochemistry?
forcModel = 'sinmod'; % What forcing?

% define forcing name (i.e. forcing subfolder)
forcName = 'FramStrait_dep0-100';

% define a run name
runName = 'FramStrait_dep0-100';

% define experiment name
expName = 'AWI-Hausgarten';

% select years
years = 2017:2018;

% set input directories and files
forcDir = fullfile('DATA', 'particle_trajectories', forcModel, forcName);
forcDummy = 'particles_MODEL_DOMAIN_EXPERIMENT-YEAR_YEAR_t*iSub03.mat';

% % set output directory
% outDir = fullfile('OUTPUT', bioModel, 'FORCSTR', runName, expName);
% outDir = strrep(outDir, 'FORCSTR', sprintf('FRC_%s', forcModel));

% list of biological forcing variables (used for initalization)
switch lower(forcModel)
    case 'biomas' % for BIOMAS forcing
        bioForcing = {'NO3', 'Si', 'PS', 'PL', 'ZS', 'ZL', 'ZP'};
    case 'sinmod' % for SINMOD forcing
        bioForcing = {'NO3', 'PS', 'PL'};
end

% select trajectores - choose all trajectories, they're filtered later
%  - selection options
%     1. use all trajectories by setting useTraj = [];
%     2. manually define trajectory index vector, e.g.: iTraj0 = 1:100:7500

useTraj = [];
% useTraj = 1:25:5000;


% loop over selected forcing files
for iy = 1:length(years)
    % load forcing
    forcFile = replace(forcDummy, {'MODEL', 'DOMAIN', 'EXPERIMENT', 'YEAR'}, ...
                                  {forcModel, forcName, expName, num2str(years(iy))});
    forcList = dir(fullfile(forcDir, forcFile));
    forcFile = forcList.name;
    try
        tStart = tic;
        load(fullfile(forcDir, forcFile));
        forcing = PDat;
        grid = runDat.grid;
        clear PDat runDat relDat stepsDat
    catch
        error('Forcing file does not exist:\n ==> %s.', fullfile(forcDir, forcFile));
    end
    
    % filter trajectories with geographical jumps
    forcing = filterTrajectoryJumps(forcing);

    % match forcing longitudes with model and mapping coordinate system
    iWest = forcing.x >= 180.0;
    forcing.x(iWest) = forcing.x(iWest) - 360.0;
    clear iWest

    % set model forcing
    nTraj = size(forcing.x,2);
    if  isempty(useTraj)
        % no specific trajectories seleced => use all
        iTraj = 1:nTraj;
    else
        iSel = false(1,nTraj);
        if  sum(useTraj>0)>0
            % used selected trajectory indices
            iSel(useTraj(useTraj>0 & useTraj<=nTraj)) = true;
        end
        iTraj = find(iSel); % convert logical mask into index vector
    end

    % filter forcing data
    y_index = ['y' num2str(years(iy))];
    F.(y_index).iTraj = iTraj;
    F.(y_index).t     = forcing.t;
    F.(y_index).z     = grid.z;
    F.(y_index).x     = forcing.x(:,iTraj);
    F.(y_index).y     = forcing.y(:,iTraj);
    F.(y_index).H     = -forcing.H(:,iTraj);
    F.(y_index).ice   = forcing.ice(:,iTraj);
    F.(y_index).iceh  = forcing.iceh(:,iTraj);
    F.(y_index).swrad = forcing.swrad(:,iTraj);
    F.(y_index).T     = forcing.profiles.temp(:,iTraj,:);
    F.(y_index).kv    = forcing.profiles.Ks(:,iTraj,:) * 0.0001; % cm2/s -> m2/s
    % The 'bioForcing' variables returned from SINMOD can be used to help
    % initialise the model
    for j = 1:length(bioForcing)
        F.(y_index).(sprintf('%sic', bioForcing{j})) = ... 
            forcing.profiles.(bioForcing{j})(:,iTraj,:);
    end
    
    clear forcing
    % Permute trajectory to last dimension
    F.(y_index).T = permute(F.(y_index).T , [1, 3, 2]);
    F.(y_index).kv = permute(F.(y_index).kv , [1, 3, 2]);        
    for j = 1:length(bioForcing)
        F.(y_index).(sprintf('%sic', bioForcing{j})) = ...
            permute(F.(y_index).(sprintf('%sic', bioForcing{j})), [1, 3, 2]);
    end
end



%% Model set-up

% Model input is stored in 4 structs, Forc, FixedParams, Params, and Data, 
% containing particle-trajectory forcing data, constant parameters,
% variable parameters that may be numerically tuned, and observed data used
% to tune parameters.

% Load/prepare in-situ fitting data
obsDir = fullfile('DATA', 'AWI_Hausgarten');
Data = prepareFittingData(obsDir, 'plotRawSizeSpectra', false, ...
    'plotSizeClassBins', false, 'plotRawSpectraAndBins', true); % tailor this function to specific data sets

% Choose model dimensions and initial parameter values. The modelled cell
% size ranges are chosen to correspond with the size class intervals
% selected in prepareFittingData.m
[FixedParams, Params] = initialiseParameters(F, Data);
Params0 = Params;

% Parameter values can be changed using name-value pair arguments in 
% updateParameters.m, eg,
% Params = updateParameters(Params, FixedParams, 'pmax_a', 30, 'pmax_b', -0.55, 'Gmax', 3);

% Interpolate forcing data over modelled depth layers
Forc = prepareForcing(F,FixedParams);
clear F

% Remove fitting-data samples from below the maximum modelled depth
Data = omitDeepSamples(Data, FixedParams);

% Choose years to model from the available data, return filtered data
[Data, Forc, FixedParams] = selectYears(Data, Forc, FixedParams);

% Choose set of forcing data trajectories near to sampling sites. 
maxDist = 25; % choose trajectories from within maxDist km radius of sampling sites
numTraj = 10; % number of trajectories per sample event (duplicates may be 
              % required for some events, depending on maxDist)
[Forc, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, numTraj);

% Omit data from any sampling events not matched with a set of trajectories
[Data, eventTraj] = omitUnmatchedEvents(Data, eventTraj, Forc);

% Find time of latest sample along each trajectory so we can reduce 
% integration lengths during parameter optimisation
Forc = latestSampleTime(Forc, Data);

% Group sampling events by origin of particles: Arctic or Atlantic.
[Forc, Data.scalar] = particleOrigin(Forc, Data.scalar, ...
    'trajectoryPlot', true, 'dendrogram', true); pause(0.25)

% Standardise the fitting data using linear mixed models to adjust for
% variability due to depth and sampling event.
Data = standardiseFittingData(Data, 'plotScaledPON', true, 'plotScaledPOC', true, ...
    'plotScaledchl_a', true, 'plotScaledN', true, 'plotScaledN_at_size', true, ...
    'plotAllData', true);
% These linear mixed models do a good job standaridsing the data. However,
% as they do not actually represent the real processes generating depth-
% and event-related discrepancies in measured values, there are a couple of
% outlying points that will cause statistical issues within the cost
% function. There are only 3 of these outliers that I'm going to remove
% (those points beyond 3 standard deviations of the mean, see plots of N
% and Chl). Remove these outliers then standardise again...
% This could be improved (perhaps fewer outliers) by altering the linear
% mixed models to account for changing curvature of measured variables
% near the surface water, as the log transform cannot capture that shape.
ind = Data.scalar.scaled_Value > 3 | Data.scalar.scaled_Value < -3; % index the outliers
fields = fieldnames(Data.scalar);
for i = 1:length(fields)
    if size(Data.scalar.(fields{i}),1) == Data.scalar.nSamples
        Data.scalar.(fields{i}) = Data.scalar.(fields{i})(~ind);
    end
    if i == length(fields)
        Data.scalar.nSamples = Data.scalar.nSamples - sum(ind);
    end
end

Data = standardiseFittingData(Data, 'plotScaledPON', false, 'plotScaledPOC', false, ...
    'plotScaledchl_a', false, 'plotScaledN', false, 'plotScaledN_at_size', false, ...
    'plotAllData', false);

% Smooth the data using polynomials to represent data shape. Fitting to the
% data shape instead of each data point should reduce effect of noise
% in the data and make the cost function more robust
Data.scalar.maxDegree = 3; % Polynomial degree = 3 should be OK if the
                           % standardised data have a monomodal distribution
% Data.size.maxDegree = 2; % Quadratic for the sorted size data seems reasonable, anything more looks like overfitting
Data.size.maxDegree = 1; % Linear for the scaled then sorted size data looks OK
Data = smoothData_polynomials(Data, Data.scalar.maxDegree, Data.size.maxDegree, ...
    'plotPolyPON', true, 'plotPolyPOC', true, 'plotPolyN', true, ...
    'plotPolyChl', true, 'plotPolySize', true);

close all


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% State variable initial condition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
% 
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [depth]
%                      organic matter      [type, depth, nutrient]
%

v0 = initialiseVariables(FixedParams, Params, Forc); % v0 stores initial input vectors


%% Integrate

poolObj = gcp; % integrations are parallelised over trajectories
numcores = poolObj.NumWorkers;

% Integrator is called separately for each time step -- defined by forcing 
% data (daily) intervals
odeMaxTime = FixedParams.dt_max; % max timestep (days)
odeInitTime = 0.5 * odeMaxTime; % initial integration timestep (solver will automatically reduce this if required)

% Choose solver
integratorChoices = {'ode45', 'ode23', 'ode113'};
odeIntegrator = integratorChoices{2}; % ode23 seems to be the most robust. ode45 occasionally produced NaNs... maybe ode45 is less robust due to stiffness in the model equations...

% Set solver options: positive definite, tolerances, initial & max time steps
odeOptions=odeset('AbsTol',1e-6,'RelTol',1e-4,...
  'InitialStep',odeInitTime,'MaxStep',odeMaxTime, ...
  'NonNegative', ones(FixedParams.nEquations,1));

% Run the model
tic
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);
runtime = toc;
disp([num2str(runtime) ' seconds to integrate ' num2str(Forc.nTraj) ' trajectories'])
runtime = runtime / Forc.nTraj * numcores; % average integration time of single trajectory on single processor

% out
% auxVars



%% Parameter tuning

clear out auxVars OUT AUXVARS AUXVARS_2d

poolObj = gcp; % start parallel pool
numcores = poolObj.NumWorkers;

integratorChoices = {'ode45', 'ode23', 'ode113'};
odeIntegrator = integratorChoices{2}; % select ODE integrator

optimiserChoices = {'ga','muga'};
optimiser = optimiserChoices{1}; % select optimising algorithm (so far only ga is available)

costFunctionChoices = {'LSS','polyLikelihood'};
costFunctionType = costFunctionChoices{2}; % select cost function
FixedParams.costFunction = costFunctionType;

% Choose tuning parameters from Params.scalars and Params.sizeDependent.
parnames = {'rDON', 'rPON', 'rPOC', 'beta2', 'beta3', 'aP', 'Gmax', 'k_G', ...
    'pmax_a', 'pmax_b', 'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', ...
    'Qmax_delQ_b', 'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b'};
npars = length(parnames);
par0 = nan(1,npars); lb = nan(1,npars); ub = nan(1,npars);
for i = 1:npars
    par0(i) = Params.(parnames{i});
    lb(i) = Params.lowerBound.(parnames{i});
    ub(i) = Params.upperBound.(parnames{i});
end
FixedParams.tunePars = parnames;
FixedParams.tunePars_lb = lb;
FixedParams.tunePars_ub = ub;

Forc.integrateFullTrajectory = false;

% Test the cost function
tic
[cost, costComponents, modData, out, auxVars] = costFunction(par0, ...
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, ... 
    'selectFunction', costFunctionType);
runtime = toc;
disp([num2str(runtime) ' seconds to integrate ' num2str(Forc.nTraj) ' trajectories'])
runtime = runtime / Forc.nTraj * numcores; % average integration time of single trajectory on single processor
disp(costComponents)

% Tune parameters
popSize = 100; % number of separate parameter sets in population of optimising algorithm
niter = 10; % algorithm iterations

if exist('runtime', 'var'), approxOptimisationTime = ... 
        approxRunTime(runtime, numcores, Forc.nTraj, popSize, niter+1, 'Display', true); end

gaOptions = gaoptimset('PopInitRange', [lb;ub], 'PopulationSize', popSize, ...
    'Generations', niter, 'Display', 'iter', ...
    'OutputFcns', @gaStoreHistory, 'PlotFcns', @gaplotbestf);

% Call optimiser
clear cost costComponents modData out auxVars optPars fval
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = ga(@(x) costFunction( ...
    x, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, ...
    'selectFunction', costFunctionType), npars, [], [], [], [], lb, ub, [], gaOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now'))


% Store results in gaOutput
stoppedEarly = ~exist('fval', 'var');
gaOutput = storeGAoutput(gapopulationhistory, gacosthistory, gaOptions, ... 
    FixedParams, 'stoppedEarly', stoppedEarly);

% Save/load output
tag = '_1';
fileName = ['results/fittedParameters_' FixedParams.costFunction tag];

saveParams = true;
switch saveParams
    case true
        saveOptimisationRun(fileName, gaOutput, Data, Forc, FixedParams, Params, v0);
end

loadParams = true;
switch loadParams
    case true
        [~, gaOutput, parnames, optPar, lb, ub, Data, Forc, FixedParams, ... 
            Params, v0] = loadOptimisationRun(fileName);
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% We may continue optimising parameters stored in gaOutput by initialising
% the algorithm with a stored population.
populationOld = gaOutput.populationHistory(:,:,end);
scoresOld = gaOutput.scoreHistory(:,end);
gaOptions = gaOutput.gaOptions;
niter = 10;

gaOptions = gaoptimset(gaOptions, 'Generations', niter, ... 
    'InitialPopulation', populationOld, 'InitialScores', scoresOld);

tic; disp('.. started at'); disp(datetime('now'))
clear optPars fval
[optPar, fval, exitflag, output, population, scores] = ga(@(x) costFunction( ...
    x, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, ... 
    'selectFunction', costFunctionType), npars, [], [], [], [], lb, ub, [], gaOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now'))

% Store results in gaOutput
stoppedEarly = ~exist('fval', 'var');
gaOutput = storeGAoutput(gapopulationhistory, gacosthistory, gaOptions, ... 
    FixedParams, 'stoppedEarly', stoppedEarly);


% Save/load output
tag = '_1';
fileName = ['results/fittedParameters_' FixedParams.costFunction tag];

saveParams = true;
switch saveParams
    case true
        saveOptimisationRun(fileName, gaOutput, Data, Forc, FixedParams, Params, v0);
end

loadParams = true;
switch loadParams
    case true
        [~, gaOutput, parnames, optPar, lb, ub, Data, Forc, FixedParams, ... 
            Params, v0] = loadOptimisationRun(fileName);
end



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Generate output using fitted parameters
Params = updateParameters(Params, FixedParams, optPar);

Forc.integrateFullTrajectory = true;

clear cost costComponents modData out auxVars
[cost,costComponents,modData,out,auxVars] = costFunction(optPar, ...
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, ... 
    'selectFunction', costFunctionType);


%% Plots

%~~~~~~~~~~~~~~~~~~
% Model fit to data
%~~~~~~~~~~~~~~~~~~

close all

save = false;
% save = true;
folder = 'results/plots/';

% Summary plots displaying model fit to data
pltChl = outputPlot('outputVsData_summaryPlots', 'chl_a', Data, modData, FixedParams); pause(0.5)
pltN = outputPlot('outputVsData_summaryPlots', 'N', Data, modData, FixedParams); pause(0.5)
pltPON = outputPlot('outputVsData_summaryPlots', 'PON', Data, modData, FixedParams); pause(0.5)
pltPOC = outputPlot('outputVsData_summaryPlots', 'POC', Data, modData, FixedParams); pause(0.5)
pltSize = outputPlot('outputVsData_summaryPlots', 'N_at_size', Data, modData, FixedParams); pause(0.5)

switch save
    case true
        filename = 'fitToData_chl.png';
        print(pltChl, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'fitToData_DIN.png';
        print(pltN, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'fitToData_PON.png';
        print(pltPON, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'fitToData_POC.png';
        print(pltPOC, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'fitToData_sizeSpectra.png';
        print(pltSize, fullfile(folder, filename), '-r300', '-dpng');
end


%~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot single trajectories
%~~~~~~~~~~~~~~~~~~~~~~~~~

% Save plots?
save = false;
% save = true;
folder = 'results/plots/';

close all

% Choose trajectory
% k = 1;
% Or, first, filter by sampling event 
ie = 20; % sampling event
kk = find(Data.scalar.EventTraj(ie,:)); % all trajectories for selected event
k = kk(5);

% Depth-time contour plots
pltForc = outputPlot('contour_DepthTime','forcing',k,out,FixedParams,Forc,auxVars,'linear');
pltIN = outputPlot('contour_DepthTime','inorganicNutrient',k,out,FixedParams,Forc,auxVars,'linear');
pltOM = outputPlot('contour_DepthTime','DOM_POM',k,out,FixedParams,Forc,auxVars,'linear');
pltPN = outputPlot('contour_DepthTime','phytoplankton_N',k,out,FixedParams,Forc,auxVars,'linear');
pltChl = outputPlot('contour_DepthTime','phytoplankton_Chl',k,out,FixedParams,Forc,auxVars,'linear');
pltPC = outputPlot('contour_DepthTime','phytoplankton_C',k,out,FixedParams,Forc,auxVars,'linear');
pltPNC = outputPlot('contour_DepthTime','phytoplankton_N_C',k,out,FixedParams,Forc,auxVars,'linear');
pltPChlN = outputPlot('contour_DepthTime','phytoplankton_Chl_N',k,out,FixedParams,Forc,auxVars,'linear');
pltZ = outputPlot('contour_DepthTime','zooplankton',k,out,FixedParams,Forc,auxVars,'linear');

switch save
    case true
        filename = 'forcing_data.png';
        print(pltForc, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'DIN.png';
        print(pltIN, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'OM.png';
        print(pltOM, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'phytoplankton_N.png';
        print(pltPN, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'phytoplankton_Chl.png';
        print(pltChl, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'phytoplankton_C.png';
        print(pltPC, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'phytoplankton_N_C_ratio.png';
        print(pltPNC, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'phytoplankton_Chl_N_ratio.png';
        print(pltPChlN, fullfile(folder, filename), '-r300', '-dpng');
        filename = 'zooplankton.png';
        print(pltZ, fullfile(folder, filename), '-r300', '-dpng');
end



%~~~~~~~~~~~~~~~
% Time evolution
%~~~~~~~~~~~~~~~

close all

% Choose event
ie = 20;
if ~ismember(ie, 1:Data.scalar.nEvents), warning(['Choose event number within range (1, ' num2str(Data.scalar.nEvents) ')']); end
% trajectory indices
kk = find(Data.scalar.EventTraj(ie,:));

outputPlot('trajectoryPolygon_TimeSeries','forcing',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','DIN',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','DOM_POM',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoplankton_C',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoplanktonStacked',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoZooPlanktonStacked',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('barplot_TimeSeries','phytoZooPlankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot groups of trajectories corresponding to each event
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

close all

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


