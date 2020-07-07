%% Size-structured 1D monod model with a single predator class

clear
clc
close all

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
% useTraj = 1:100:5000;


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

% Model input is stored in 3 structs, 'Forc', 'FixedParams' and 'Params',
% containing particle-trajectory forcing data, constant parameters, and
% variable parameters that may be numerically tuned. A 4th struct, 'Data',
% containing measurements is required if parameters are numerically tuned.

% Choose model dimensions and initial parameter values
[FixedParams, Params] = initialiseParameters(F);
Params0 = Params;

% Parameter values can be changed using name-value pair arguments in 
% updateParameters.m, eg,
% Params = updateParameters(Params, FixedParams, 'pmax_a', 30, 'pmax_b', -0.55, 'Gmax', 3);


%~~~~~~~~~~~~~
% Prepare data
%~~~~~~~~~~~~~

% Interpolate forcing data over chosen depth layers
Forc0 = prepareForcing(F,FixedParams);
% clear F

% Load and prepare fitting data
obsDir = fullfile('DATA', 'AWI_Hausgarten');
Data = prepareFittingData(obsDir, FixedParams); % tailor this function to specific data sets

% Use forcing data from years corresponding to samples
iy = ismember(Forc0.years, unique(Data.Year));
disp(['use particle trajectories from ' num2str(Forc0.years(iy))])
% Filter out unused years
[Y,~] = datevec(Forc0.t);
ind = any(Y == Forc0.years(iy));
ntraj = length(ind);
ntraj_new = sum(ind);
fields = fieldnames(Forc0);
for i = 1:length(fields)
    x = Forc0.(fields{i});
    s = size(x);
    if any(s == ntraj)
        nd = length(s);
        ind_ = repmat(reshape(ind, [ones(1,nd-1) ntraj]), [s(1:end-1) 1]);
        Forc.(fields{i}) = reshape(x(ind_), [s(1:end-1) ntraj_new]);
    else
        Forc.(fields{i}) = Forc0.(fields{i});
    end
end
Forc.years = Forc.years(iy);
FixedParams.years = FixedParams.years(iy);


% Filter forcing data by finding trajectories close to sampling events
maxDist = 25; % distance of particle from sample location
maxTraj = 10; % maximum number of particles per sample event
[Forc, eventTraj] = chooseTrajectories(Forc, Data, maxDist, maxTraj);

% If any sampling events have not been ascribed trajectories then either
% refilter using above code, or omit those events
if Data.nEvents > length(unique(eventTraj.event))
    Data.nEvents = length(unique(eventTraj.event));
    keep = ismember(Data.Event, eventTraj.event);
    fields = fieldnames(Data);
    for i = 1:length(fields)
        if size(Data.(fields{i}), 1) == Data.nSamples
            Data.(fields{i}) = Data.(fields{i})(keep);
        end
    end
    Data.nSamples = size(Data.Value,1);
    adj = [unique(Data.Event) (1:Data.nEvents)'];
    for i = 1:size(adj,1)
        Data.Event(Data.Event == adj(i,1)) = adj(i,2);
        eventTraj.event(eventTraj.event == adj(i,1)) = adj(i,2);
    end
end    

Data.EventTraj = false(Data.nEvents, Forc.nTraj); % match samples to trajectories
for i = 1:Data.nEvents
    ti = eventTraj.trajIndex(eventTraj.event == i);
    Data.EventTraj(i,ti) = true;
end


% Group sampling events by origin of particles: Arctic or Atlantic.
[Forc, Data] = particleOrigin(Forc, Data, ...
    'trajectoryPlot', true, 'dendrogram', true);

% Standardise the fitting data using linear mixed models to adjust for
% variability due to depth and sampling event.
Data = standardiseFittingData(Data,'plotScaledPON', true, 'plotScaledPOC', true, ...
    'plotScaledchl_a', true, 'plotScaledN', true, 'plotAllData', true);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% State variable initial condition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
% 
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth]
%                      zooplankton         [depth]
%                      organic matter      [depth]
%

v0 = initialiseVariables(FixedParams,Forc); % v0 stores initial input vectors


%% Integrate

% Solve with ode45, called separately for each time step -- defined by
% forcing data (daily) intervals
odeMaxTime = FixedParams.dt_max; % max timestep (days)
odeInitTime = 0.5 * odeMaxTime; % initial timestep (ode45 will automatically reduce this if required)

% set solver options: positive definite, tolerances, initial & max time steps
ode45options=odeset('NonNegative',[1 ones(1, FixedParams.nEquations)],...
  'AbsTol',1e-6,'RelTol',1e-4,...
  'InitialStep',odeInitTime,'MaxStep',odeMaxTime);

poolObj = gcp; % integrations are parallelised over trajectories
numcores = poolObj.NumWorkers;

% Integrate model specified by Params
tic
[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options);
runtime = toc;
disp([num2str(runtime) ' seconds to integrate ' num2str(Forc.nTraj) ' trajectories'])
runtime = runtime / Forc.nTraj * numcores; % average integration time of single trajectory on single processor

% Extract solutions from array, OUT, into a more readable struct, out.
% Same for extra outputs...
nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nOM = FixedParams.nOM;
nTraj = Forc.nTraj;

out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM nz nt nTraj]);

if sum(nExtra) > 0
    for k = 1:nExtra(1)
        auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
    end
    for k = 1:nExtra(2)
        auxVars.(namesExtra{k+nExtra(1)}) = squeeze(AUXVARS_2d(:,:,k,:,:));
    end
end


%% Parameter tuning

clear OUT out AUXVARS AUXVARS_2d auxVars

poolObj = gcp; % integrations are parallelised over trajectories
numcores = poolObj.NumWorkers;

% Choose tuning parameters from Params.scalars and Params.sizeDependent
parnames = {'pmax_a','pmax_b','aP', ... 
    'Gmax','k_G', ... 
    'aN_over_Qmin_a','aN_over_Qmin_b'};
npars = length(parnames);
FixedParams.tunePars = parnames;
lb = nan(1,npars); ub = nan(1,npars);
for i = 1:npars
    lb(i) = Params.lowerBound.(parnames{i});
    ub(i) = Params.upperBound.(parnames{i});
end

approxRunTime = @(runtime, ncores, ntraj, popsize, niter) ... % assuming integrating over a single trajectory takes 1 sec on a single processor
    round(runtime * ntraj * (niter+1) * popsize / ncores / 60 / 60, 2, 'significant');

popSize = 50;
niter = 10;

approxOptimisationTime = approxRunTime(runtime, numcores, Forc.nTraj, popSize, niter);
disp(['predicted optimsation run-time ~ ' num2str(approxOptimisationTime) ' hrs'])

gaOptions = optimoptions('ga','Display','iter','PlotFcn',@gaplotbestf, ...
    'MaxGenerations',niter,'PopulationSize',popSize);

% tune parameters
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = ... 
    ga(@(x)costFun(x, FixedParams, Params, Forc, Data, v0, ode45options), ...
    npars, [], [], [], [], lb, ub, [], gaOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now'))

% % continue parameter-tuning from where the algorithm stopped...
% niter = 10;
% gaOptions = optimoptions('ga','Display','iter','PlotFcn',@gaplotbestf, ...
%     'MaxGenerations',niter,'InitialPopulation',population);
% tic; disp('.. started at'); disp(datetime('now'))
% [optPar, fval, exitflag, output, population, scores] = ... 
%     ga(@(x)costFun(x, FixedParams, Params, Forc, Data, v0, ode45options), ...
%     npars, [], [], [], [], lb, ub, [], gaOptions);
% optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now'))


% generate output using fitted parameters
Params = updateParameters(Params, FixedParams, optPar);

[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options);

% Extract solutions
nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nOM = FixedParams.nOM;
nTraj = Forc.nTraj;
out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM nz nt nTraj]);
if sum(nExtra) > 0
    for k = 1:nExtra(1)
        auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
    end
    for k = 1:nExtra(2)
        auxVars.(namesExtra{k+nExtra(1)}) = squeeze(AUXVARS_2d(:,:,k,:,:));
    end
end



%% Plot some output...

%~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot single trajectories
%~~~~~~~~~~~~~~~~~~~~~~~~~

% Save plots?
save = false;
% save = true;
folder = 'OUTPUT/plots/';

close all

% Choose trajectory
% k = 1;
% Or, first, filter by sampling event 
ie = 10; % sampling event
kk = find(Data.EventTraj(ie,:)); % all trajectories for selected event
k = kk(1);

% Depth-time contour plots

% forcing data
outputPlot('contour_DepthTime','forcing',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'forcing_data.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% inorganic nutrient
outputPlot('contour_DepthTime','inorganicNutrient',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'dissolved_inorganic_nutrient.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% organic nutrient
outputPlot('contour_DepthTime','DOM_POM',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'organic_nutrient.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton
outputPlot('contour_DepthTime','phytoplankton',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton carbon
outputPlot('contour_DepthTime','phytoplankton_C',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton N/C ratio
outputPlot('contour_DepthTime','phytoplankton_N_C',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% zooplankton
outputPlot('contour_DepthTime','zooplankton',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'zooplankton.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot groups of trajectories corresponding to each event
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

close all

% Choose event
ie = 2;
if ~ismember(ie, 1:Data.nEvents), warning(['Choose event number within range (1, ' num2str(Data.nEvents) ')']); end
% trajectory indices
kk = find(Data.EventTraj(ie,:));

outputPlot('trajectoryLine_LatLong','direction',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','forcing',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','inorganicNutrient',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','DOM_POM',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','phytoplankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','zooplankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);



%~~~~~~~~~~~~~~~
% Time evolution
%~~~~~~~~~~~~~~~

close all

% Choose event
ie = 1;
if ~ismember(ie, 1:Data.nEvents), warning(['Choose event number within range (1, ' num2str(Data.nEvents) ')']); end
% trajectory indices
kk = find(Data.EventTraj(ie,:));

outputPlot('trajectoryPolygon_TimeSeries','forcing',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','DOM_POM',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoplankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoplanktonStacked',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoZooPlanktonStacked',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);

outputPlot('barplot_TimeSeries','phytoZooPlankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);


