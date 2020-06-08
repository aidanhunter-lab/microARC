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

% Parameter values can be changed using name-value pairs in
% updateParameters.m, eg,
% Params = updateParameters(Params, FixedParams, 'pmax_a', 30, 'pmax_b', -0.55, 'Gmax', 3);


%~~~~~~~~~~~~~
% Prepare data
%~~~~~~~~~~~~~

% Interpolate forcing data over chosen depth layers
F = prepareForcing(F,FixedParams);

% load fitting data
obsDir = fullfile('DATA', 'AWI_Hausgarten');
obsFile = 'AWI-Hausgarten.xlsx';
dat = readtable(fullfile(obsDir,obsFile));
% rename variables
vn = dat.Properties.VariableNames;
vn{strcmp(vn, 'PangaeaEventLabel')} = 'eventLabel';
vn{strcmp(vn, 'lat_degN')} = 'lat';
vn{strcmp(vn, 'long_degE')} = 'lon';
vn{strcmp(vn, 'Sample')} = 'sample';
vn{strcmp(vn, 'depth_m')} = 'depth';
vn{strcmp(vn, 'StationName')} = 'station';
dat.Properties.VariableNames = vn;
% use years with samples and forcing data
clear year
iForc = ismember(year(dat.SamplingDate), FixedParams.years);
dat = dat(iForc,:);
dat.year = year(dat.SamplingDate);
dat.month = month(dat.SamplingDate);
dat.day = day(dat.SamplingDate);
dat.yearday = floor(yearday(datenum(dat.SamplingDate)));
eventLabel = unique(dat.eventLabel, 'stable');
event = (1:length(eventLabel))';
events = table(event, eventLabel);
dat = join(dat, events);
% remove duplicate event numbers corresponding to samples at different depths
dat0 = dat;
dat.depth = [];
dat.sample = [];
dat = unique(dat, 'stable');
head(dat)

% Filter forcing data by finding trajectories close to sampling events
maxDist = 25; % distance of particle from sample location
maxTraj = 10; % maximum number of particles per sample event
[Forc, eventTraj] = chooseTrajectories(F, dat, maxDist, maxTraj);
% head(eventTraj) % list all trajectories selected for each sampling event

% Store fitting-data in separate struct
events = unique(eventTraj.event);
Data.nEvent = length(events);
Data.label = dat.eventLabel(ismember(dat.event,events));
Data.t = datenum(dat.SamplingDate(ismember(dat.event,events)));
Data.eventTraj = false(Data.nEvent, Forc.nTraj);
for i = 1:Data.nEvent
    ti = eventTraj.trajIndex(eventTraj.event == events(i));
    Data.eventTraj(i,ti) = true;
end


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


%~~~~~~~~~~~~~~~~~~~~~~~
% Optional extra outputs
%~~~~~~~~~~~~~~~~~~~~~~~

% Can return extra output along with the integrated state variables. Four
% options (so far) for returnExtras.
% returnExtras = 'none';
returnExtras = 'auxiliary'; % return non-state variables
% returnExtras = 'rates';   % return rates of change of state variables
% returnExtras = 'auxiliaryAndRates';
FixedParams.returnExtras = returnExtras;


%% Integrate

% Solve with ode45, called separately for each time step -- defined by
% forcing data (daily) intervals
odeInitTime = 12/24; % initial timestep = 12 hrs. (ode45 will automatically reduce this if required)
odeMaxTime = 1;      % max timestep = 1 day

% set solver options: positive definite, tolerances, initial & max time steps
ode45options=odeset('NonNegative',[1 ones(1, FixedParams.nEquations)],...
  'AbsTol',1e-6,'RelTol',1e-4,...
  'InitialStep',odeInitTime,'MaxStep',odeMaxTime);

poolObj = gcp; % integrations are parallelised over trajectories


% integrating function
tic
[OUT, AUXVARS, AUXVARS_2d, RATES, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options);
toc


% %---------------------------------
% % Parameters to tune
% parnames = {'A', 'm', 'aP', 'Gmax', 'k_G', 'Lambda', 'lambda_max', 'wk', ...
%     'rPOM', 'rDOM', 'Qmin_a', 'Qmin_b', 'Qmax_over_delQ_a', ... 
%     'Qmax_over_delQ_a', 'Vmax_over_Qmin_a', 'Vmax_over_Qmin_b', ... 
%     'aN_over_Qmin_a', 'aN_over_Qmin_b', 'pmax_a', 'pmax_b'};
% FixedParams.tunePars = parnames;
% npars = length(parnames);
% pars = nan(1, npars);
% for i = 1:npars
%     pars(i) = Params.(parnames{i});
% end
% 
% % tune parameters
% ga(@(x)costFun(x, FixedParams, Params, Forc, Data, v0, ode45options), npars)
% 
% %---------------------------------


% Extract solutions from array, OUT, into a more readable struct, out.
% Same for extra outputs...
nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nTraj = Forc.nTraj;

out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [1 nz nt nTraj]);

if ~strcmp(returnExtras, 'none')
    switch returnExtras
        case 'auxiliary'
            for k = 1:nExtra(1)
                auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
            end
            for k = 1:nExtra(2)
                auxVars.(namesExtra{nExtra(1) + k}) = squeeze(AUXVARS_2d(:,:,k,:,:));
            end
        case 'auxiliaryAndRates'
            for k = 1:nExtra(1)
                auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
            end
            for k = 1:nExtra(2)
                auxVars.(namesExtra{nExtra(1) + k}) = squeeze(AUXVARS_2d(:,:,k,:,:));
            end
            rates.N = reshape(RATES(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
            rates.P = reshape(RATES(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
            rates.Z = reshape(RATES(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
            rates.OM = reshape(RATES(FixedParams.OM_index,:,:), [1 nz nt nTraj]);
        case 'rates'
            rates.N = reshape(RATES(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
            rates.P = reshape(RATES(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
            rates.Z = reshape(RATES(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
            rates.OM = reshape(RATES(FixedParams.OM_index,:,:), [1 nz nt nTraj]);
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
k = 1;
% Or, first, filter by year or sampling event 
iy = 1;  % year index
ie = 10; % sampling event

if FixedParams.years(iy) ~= dat.year(dat.event == ie)
    kk = nan;
    warning('Selected sampling event must occur during selected year')
else
    kk = eventTraj.trajIndex(eventTraj.event == ie); % all trajectories for selected event
end
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
outputPlot('contour_DepthTime','organicNutrient',k,out,FixedParams,Forc,auxVars,'linear');
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
ie = 8;
if ~ismember(ie, eventTraj.event), warning(['Choose event number within range (' num2str(min(eventTraj.event)) ', ' num2str(max(eventTraj.event)) ')']); end
% trajectory indices
kk = eventTraj.trajIndex(eventTraj.event == ie);

outputPlot('trajectoryLine_LatLong','direction',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryLine_LatLong','forcing',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryLine_LatLong','inorganicNutrient',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryLine_LatLong','organicNutrient',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryLine_LatLong','phytoplankton',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryLine_LatLong','zooplankton',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);



%~~~~~~~~~~~~~~~
% Time evolution
%~~~~~~~~~~~~~~~

close all

% Choose event
ie = 1;
if ~ismember(ie, eventTraj.event), warning(['Choose event number within range (' num2str(min(eventTraj.event)) ', ' num2str(max(eventTraj.event)) ')']); end
% trajectory indices
kk = eventTraj.trajIndex(eventTraj.event == ie);

outputPlot('trajectoryPolygon_TimeSeries','forcing',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryPolygon_TimeSeries','nutrient',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoplankton',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoplanktonStacked',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);
outputPlot('trajectoryPolygon_TimeSeries','phytoZooPlanktonStacked',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);

outputPlot('barplot_TimeSeries','phytoZooPlankton',ie,kk,out,FixedParams,Forc,auxVars,dat,0.1);






