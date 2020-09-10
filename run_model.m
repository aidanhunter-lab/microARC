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

% Load/prepare fitting data
obsDir = fullfile('DATA', 'AWI_Hausgarten');
Data = prepareFittingData(obsDir, 'plotRawSizeSpectra', false, ...
    'plotSizeClassBins', false, 'plotRawSpectraAndBins', true); % tailor this function to specific data sets

% Choose model dimensions and initial parameter values. The modelled cell
% sizes are chosen to correspond with the size class intervals selected in
% prepareFittingData.m
[FixedParams, Params] = initialiseParameters(F, Data);
Params0 = Params;

% Parameter values can be changed using name-value pair arguments in 
% updateParameters.m, eg,
% Params = updateParameters(Params, FixedParams, 'pmax_a', 30, 'pmax_b', -0.55, 'Gmax', 3);
% v= {'pmax_a', 30, 'pmax_b', -0.55, 'Gmax', 3}

% Remove fitting-data samples from below the maximum modelled depth
ind = Data.scalar.Depth < -min(FixedParams.z);
if ~all(ind)    
    fields = fieldnames(Data.scalar);
    for i = 1:length(fields)
        if size(Data.scalar.(fields{i}), 1) > 1
            Data.scalar.(fields{i}) = Data.scalar.(fields{i})(ind);
        end
    end
    Data.scalar.nSamples = length(Data.scalar.Year);
    Data.scalar.nEvents = length(unique(Data.scalar.Event));
end
ind = Data.size.DepthMin < -min(FixedParams.z);
if ~all(ind)
    fields = fieldnames(Data.size);
    for i = 1:lenth(fields)        
        if size(Data.size.(fields{i}), 1) > 1        
            Data.size.(fields{i}) = Data.size.(fields{i})(ind);
        end
    end
    Data.size.nSamples = length(Data.size.Year);
end


% Interpolate forcing data over modelled depth layers
Forc0 = prepareForcing(F,FixedParams);
clear F


% Use forcing data from years corresponding to samples
nut_years = unique(Data.scalar.Year(strcmp(Data.scalar.Variable, 'N')));
OM_years = unique(Data.scalar.Year(strcmp(Data.scalar.Variable, 'PON') | ...
    strcmp(Data.scalar.Variable, 'POC') | ... 
    strcmp(Data.scalar.Variable, 'chl_a')));
size_years = unique(Data.size.Year(strcmp(Data.size.Variable, 'N_at_size')));
use_years = size_years(ismember(size_years, unique([nut_years OM_years])));
iy = ismember(Forc0.years, use_years);
disp(['use particle trajectories from ' num2str(Forc0.years(iy))])

% Filter out unused years of forcing and/or fitting data
[Y,~] = datevec(Forc0.t);
ind = any(ismember(Y, Forc0.years(iy)));
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

% With these data available it makes sense to use 2017 data for everything,
% except the organic matter for which we only have 2016 data.
ind = ismember(Data.size.Year, Forc.years);
fields = fieldnames(Data.size);
for i = 1:length(fields)
    if size(Data.size.(fields{i}), 1) == Data.size.nSamples
       Data.size.(fields{i}) = Data.size.(fields{i})(ind); 
    end
end
Data.size.nSamples = length(Data.size.Year);
% Relabel the event numbers
eventLabs = table((1:Data.scalar.nEvents)', unique(Data.scalar.Event));
eventLabs.Properties.VariableNames = {'newEvent', 'Event'};
tmp = table(Data.scalar.Event); tmp.Properties.VariableNames = {'Event'};
tmp = join(tmp, eventLabs);
Data.scalar.Event = tmp.newEvent;
clear eventLabs tmp

% Filter forcing data by finding trajectories close to sampling events
maxDist = 25; % distance of particle from sample location
maxTraj = 10; % maximum number of particles per sample event
[Forc, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, maxTraj);


% If any sampling events have not been ascribed trajectories then either
% refilter using above code, or omit those events
if Data.scalar.nEvents > length(unique(eventTraj.event))
    Data.scalar.nEvents = length(unique(eventTraj.event));
    keep = ismember(Data.scalar.Event, eventTraj.event);
    fields = fieldnames(Data.scalar);
    for i = 1:length(fields)
        if size(Data.scalar.(fields{i}), 1) == Data.scalar.nSamples
            Data.scalar.(fields{i}) = Data.scalar.(fields{i})(keep);
        end
    end
    Data.scalar.nSamples = size(Data.scalar.Value,1);
    adj = [unique(Data.scalar.Event) (1:Data.scalar.nEvents)'];
    for i = 1:size(adj,1)
        Data.scalar.Event(Data.scalar.Event == adj(i,1)) = adj(i,2);
        eventTraj.event(eventTraj.event == adj(i,1)) = adj(i,2);
    end
end    

Data.scalar.EventTraj = false(Data.scalar.nEvents, Forc.nTraj); % match samples to trajectories
for i = 1:Data.scalar.nEvents
    ti = eventTraj.trajIndex(eventTraj.event == i);
    Data.scalar.EventTraj(i,ti) = true;
end


% Group sampling events by origin of particles: Arctic or Atlantic.
[Forc, Data.scalar] = particleOrigin(Forc, Data.scalar, ...
    'trajectoryPlot', true, 'dendrogram', true);

% Standardise the fitting data using linear mixed models to adjust for
% variability due to depth and sampling event.
% Data = standardiseFittingData(Data,'plotScaledPON', true, 'plotScaledPOC', true, ...
%     'plotScaledchl_a', true, 'plotScaledN', true, 'plotAllData', true);
Data = standardiseFittingData(Data,'plotScaledPON', true, 'plotScaledPOC', true, ...
    'plotScaledchl_a', true, 'plotScaledN', true, 'plotScaledN_at_size', true, ...
    'plotAllData', true);

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
v0 = log(v0); % integrate variables on log-scale to eliminate need for positivity constraint in odeOptions

%% Integrate

% Integrator is called separately for each time step -- defined by forcing 
% data (daily) intervals
odeMaxTime = FixedParams.dt_max; % max timestep (days)
odeInitTime = 0.5 * odeMaxTime; % initial timestep (ode45 will automatically reduce this if required)

% Choose solver
integratorChoices = {'ode45', 'ode23', 'ode113'};
odeIntegrator = integratorChoices{2}; % ode23 seems to be the most robust. ode45 occasionally produced NaNs... I think there is moderate stiffness in the model equations that ode45 can struggle to overcome...

% set solver options: positive definite, tolerances, initial & max time steps
odeOptions=odeset('AbsTol',1e-6,'RelTol',1e-4,...
  'InitialStep',odeInitTime,'MaxStep',odeMaxTime);

poolObj = gcp; % integrations are parallelised over trajectories
numcores = poolObj.NumWorkers;

% Integrate model specified by Params
tic
[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);
runtime = toc;
disp([num2str(runtime) ' seconds to integrate ' num2str(Forc.nTraj) ' trajectories'])
runtime = runtime / Forc.nTraj * numcores; % average integration time of single trajectory on single processor

OUT = exp(OUT); % exponentiate variables to natural scale

% Extract solutions from array, OUT, into a more readable struct, out.
% Same for extra outputs...
nt = FixedParams.nt;
nz = FixedParams.nz;
% nPP = FixedParams.nPP;
nPP_size = FixedParams.nPP_size;
nPP_nut = FixedParams.nPP_nut;
% nOM = FixedParams.nOM;
nOM_type = FixedParams.nOM_type;
nOM_nut = FixedParams.nOM_nut;
nTraj = Forc.nTraj;

out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP_size nz nPP_nut nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM_type nz nOM_nut nt nTraj]);

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

integratorChoices = {'ode45', 'ode23', 'ode113'};
odeIntegrator = integratorChoices{2};

% Choose tuning parameters from Params.scalars and Params.sizeDependent
parnames = {'pmax_a','pmax_b','aP', ... 
    'Gmax','k_G', ... 
    'aN_QC_a','aN_QC_b'};
npars = length(parnames);
FixedParams.tunePars = parnames;
lb = nan(1,npars); ub = nan(1,npars);
for i = 1:npars
    lb(i) = Params.lowerBound.(parnames{i});
    ub(i) = Params.upperBound.(parnames{i});
end

% Test the cost function
xpar = nan(1,npars);
for i = 1:npars
    xpar(i) = Params.(parnames{i});
end
[cost, costComponents] = costFun(xpar, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions);
disp(costComponents)

approxRunTime = @(runtime, ncores, ntraj, popsize, niter) ... % assuming integrating over a single trajectory takes 1 sec on a single processor
    round(runtime * ntraj * (niter+1) * popsize / ncores / 60 / 60, 2, 'significant');

popSize = 100;
niter = 2;

approxOptimisationTime = approxRunTime(runtime, numcores, Forc.nTraj, popSize, niter);
disp(['predicted optimsation run-time ~ ' num2str(approxOptimisationTime) ' hrs'])

gaOptions = optimoptions('ga', 'Display', 'iter', 'MaxGenerations', niter, ...
    'PopulationSize', popSize, ...
    'OutputFcn', @gaStoreHistory, ...
    'PlotFcn',@gaplotbestf);

% tune parameters
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = ... 
    ga(@(x)costFun(x, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions), ...
    npars, [], [], [], [], lb, ub, [], gaOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now'))


% Save results
gaOutput.parNames = parnames;
gaOutput.lowerBound = lb;
gaOutput.upperBound = ub;
gaOutput.optPar = optPar;
gaOutput.fval = fval;
gaOutput.exitflag = exitflag;
gaOutput.output = output;
gaOutput.population = population;
gaOutput.scores = scores;
gaOutput.populationHistory = gapopulationhistory;
gaOutput.scoreHistory = gacosthistory;

% Save output - TAKE CARE NOT TO OVER-WRITE!
m = matfile('results/fittedParameters', 'Writable', true);
m.gaOutput = gaOutput;
% Load output
m = matfile('results/fittedParameters', 'Writable', true);
gaOutput = m.gaOutput;
populationHistory = gaOutput.populationHistory;
scoreHistory = gaOutput.scoreHistory;
optPar = gaOutput.optPar;

% %%%%%%
% % debugging... if NaNs or complex numbers appeara in output
% 
% % extract iteration causing error in optimisation
% pH = populationHistory(:,:,end);
% sH = scoreHistory(:,end);
% 
% ind_err = isnan(sH);
% pErr = pH(ind_err,:); % parameters causing error in model or cost function
% 
% Params = updateParameters(Params, FixedParams, pErr);
% % integrate using problematic parameters
% integratorChoices = {'ode45', 'ode23', 'ode113'};
% odeIntegrator = integratorChoices{3};
% [OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
%     integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);
% 
% out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
% out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP_size nz nPP_nut nt nTraj]);
% out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
% out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM_type nz nOM_nut nt nTraj]);
% if sum(nExtra) > 0
%     for k = 1:nExtra(1), auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:)); end
%     for k = 1:nExtra(2), auxVars.(namesExtra{k+nExtra(1)}) = squeeze(AUXVARS_2d(:,:,k,:,:)); end
% end
% 
% % find trajectories producing invalid outputs
% fields = fieldnames(out);
% for i = 1:length(fields)
%     x = out.(fields{i}); xs = size(x);
%     x2d = reshape(x, [prod(xs(1:end-1)) xs(end)]);
%     nan_traj.(fields{i}) = find(any(isnan(x2d)));
%     neg_traj.(fields{i}) = find(any(x2d < 0));
% end
% 
% disp(nan_traj)
% disp(neg_traj)
% 
% % search for source of nans
% traj = 115;
% for i = 1:length(fields)
%     x = out.(fields{i}); xs = size(x);
%     x2d = reshape(x, [prod(xs(1:end-1)) xs(end)]);
%     ind = false(size(x2d));
%     ind(:,traj) = true; ind = reshape(ind, xs);
%     out_nan.(fields{i}) = reshape(x(ind), [xs(1:end-1)]);
% end
% % find timestep when nans first appear
% time = nan(1,length(fields));
% for i = 1:length(fields)
%     x = out_nan.(fields{i}); xs = size(x);
%     x2d = reshape(x, [prod(xs(1:end-1)) xs(end)]);
%     ind = isnan(x2d);
%     ind = [zeros(size(ind,1),1) diff(ind,[],2)];
%     time(i) = find(any(ind), 1);
% end
% disp(min(time))
% 
% %%%%%%



% continue parameter-tuning from where the algorithm stopped... NEED TO
% CHECK THIS ALL WORKS PROPERLY - MADE SOME CHANGES SINCE LAST USED...
% niter = 1;
% gaOptions = optimoptions('ga', 'Display', 'iter', 'MaxGenerations', niter, ...
%     'PlotFcn', @gaplotbestf, ...
%     'OutputFcn', @gaStoreHistory, ...
%     'InitialPopulation', gaOutput.population);
% tic; disp('.. started at'); disp(datetime('now'))
% [optPar, fval, exitflag, output, population, scores] = ... 
%     ga(@(x)costFun(x, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions), ...
%     npars, [], [], [], [], lb, ub, [], gaOptions);
% optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now'))


% generate output using fitted parameters
Params = updateParameters(Params, FixedParams, optPar);

[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);

OUT = exp(OUT); % exponentiate to natural scale

% Extract solutions
nt = FixedParams.nt; nz = FixedParams.nz;
nPP_size = FixedParams.nPP_size;
nPP_nut = FixedParams.nPP_nut;
nOM_type = FixedParams.nOM_type;
nOM_nut = FixedParams.nOM_nut;
nTraj = Forc.nTraj;

out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP_size nz nPP_nut nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM_type nz nOM_nut nt nTraj]);

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
% folder = '/home/aidan/Desktop/temp/outputPlots/';

close all

% Choose trajectory
% k = 1;
% Or, first, filter by sampling event 
ie = 7; % sampling event
kk = find(Data.scalar.EventTraj(ie,:)); % all trajectories for selected event
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

% phytoplankton nitrogen
outputPlot('contour_DepthTime','phytoplankton_N',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton_N.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton chlorophyll
outputPlot('contour_DepthTime','phytoplankton_Chl',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton_Chl.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton carbon
outputPlot('contour_DepthTime','phytoplankton_C',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton_C.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton N/C ratio
outputPlot('contour_DepthTime','phytoplankton_N_C',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton_N_C_ratio.png';
    figFile = fullfile(folder, filename);    
    print(fig, figFile, '-r300', '-dpng');
end

% phytoplankton Chl/N ratio
outputPlot('contour_DepthTime','phytoplankton_Chl_N',k,out,FixedParams,Forc,auxVars,'linear');
if save
    fig = gcf;
    filename = 'phytoplankton_Chl_N_ratio.png';
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



%~~~~~~~~~~~~~~~
% Time evolution
%~~~~~~~~~~~~~~~

close all

% Choose event
ie = 1;
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


