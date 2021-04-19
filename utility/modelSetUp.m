function [Forc, FixedParams, Params, Data] = modelSetUp(Directories, varargin)

% modelSetup.m is a wrapper for numerous functions used to load and
% organise all necessary (forcing and fitting) data, and to initialise
% model parameters.

% Model input is stored in 4 structs: Forc, FixedParams, Params, and Data,
% returned by modelSetup.m
% Forc        = forcing data from physical model particle trajectories.
% FixedParams = constant model parameters (not numerically tuned).
% Params      = model parameters that MAY be tuned.
% Data        = observed data included in cost function to tune parameters.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extractVarargin(varargin) % assign extra arguments to the workspace

extractStruct(Directories) % move fields of Directories to the workspace

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Select trajectories to use from physical model
% - options: 1. Use all trajectories by setting useTraj = []. Trajectories
%               are filtered later during model set-up.
%            2. Save RAM by manually defining trajectory index vector, 
%                e.g., useTraj = 1:50:5000 to extract every 50th trajectory
if ~exist('useTraj', 'var')
    useTraj = [];
end

% Specify years of available forcing data
if ~exist('years', 'var')
    years = 2017:2018;
end

% Load and extract relevant forcing data from physical model 'forcModel'.
F = loadForcing(Directories, years, useTraj);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Load/prepare in-situ fitting data.
if ~exist('plotCellConcSpectra', 'var')
   plotCellConcSpectra = false; 
end
if ~exist('plotBioVolSpectra', 'var')
   plotBioVolSpectra = false; 
end

% prepareFittingData.m must be tailored to specific data sets.
Data = prepareFittingData(obsDir, ...
    'plotCellConcSpectra', plotCellConcSpectra, ...
    'plotBioVolSpectra', plotBioVolSpectra);

% Choose number of modelled size class intervals using function
% chooseSizeClassIntervals.m
% min/max sizes to retain from the data
if ~exist('ESDmin', 'var')
    ESDmin = 1;
end
if ~exist('ESDmax', 'var')
    ESDmax = 200;
end

viewSizeClassChoices = false; % set true to view various data partitions to choose between
switch viewSizeClassChoices, case true
    % min/max number of modelled size classes
    if ~exist('nsizesMin', 'var')
        nsizesMin = 4; % min number of modelled size classes
    end
    if ~exist('nsizesMax', 'var')
        nsizesMax = 12; % max number of modelled size classes
    end
    sizeData = setSizeClassIntervals(Data.size, ...
        'datFull', Data.sizeFull, 'ESDmin', ESDmin, 'ESDmax', ESDmax, ...
        'nsizes', [], 'nsizesMin', nsizesMin, 'nsizesMax', nsizesMax, ...
        'plotSizeClassIntervals', true);
    disp(sizeData)
end

if ~exist('plotSizeClassIntervals', 'var')
    plotSizeClassIntervals = false;
end

% Choose number of modelled size classes
% setVariable('nsizes', 8)
if ~exist('nsizes', 'var')
   nsizes = 8;
end
% Use the same size classes for autotrophs and heterotrophs to ensure that
% grazing pressure is unbiased across sizes

% Initialise model parameters.
% Values can be modified in defaultParameters.m, which is called from
% within initialiseParameters.m.

% A file-path for saved parameters may be specified in setDirectories.m
% If a pre-existing parameter set is not loaded then set initialParamDir = [].
if ~exist('initialParamDir', 'var')
    initialParamDir = [];
end
[FixedParams, Params] = initialiseParameters(F, bioModel, ... 
    'ESDmin', ESDmin, 'ESDmax', ESDmax, 'nsizes', nsizes,  ...
    'parFile', initialParamDir);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameter values should not be changed by directly modifying Params
% because the size-dependent vector parameters will not update. Instead, 
% parameter values should be changed using name-value pair arguments in 
% updateParameters.m, e.g.,
% Params = updateParameters(Params, FixedParams, ... 
%     'pmax_a', 30, 'pmax_b', -0.55, 'k_G', 3);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Aggregate size spectra data according to size classes specified in
% FixedParams -- store in struct Data.size.dataBinned
[Data.size, Data.sizeFull] = setSizeClassIntervals(Data.size, ... 
    'datFull', Data.sizeFull, 'nsizes', nsizes, 'FixedParams', FixedParams, ...
    'plotSizeClassIntervals', plotSizeClassIntervals);

% Choose data types to use in cost function - I think bio-volume is the
% best choice for the size data
if ~exist('scalarFittingData', 'var')
    scalarFittingData = {'N','chl_a','PON','POC'};
end
if ~exist('vectorFittingData', 'var')
    vectorFittingData = {'BioVol'};
end

Data = selectCostFunctionData(Data, scalarFittingData, vectorFittingData);

clear sizeData

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Interpolate forcing data over modelled depth layers, and combine multiple
% years of forcing data into a single structure using prepareForcing.m.
% A few extra useful metrics are also calculated here.
Forc = prepareForcing(F,FixedParams);
clear F

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

% Run model over a single year? If false, then multiple years of forcing 
% data MAY be used depending on fitting-data availability
if ~exist('useSingleYear', 'var'), useSingleYear = true; end
% There are 2 size spectra for 2018 -- use only the Polarstern samples if true
if ~exist('useSingleSizeSpectra', 'var'), useSingleSizeSpectra = true; end

[Data, Forc, FixedParams] = selectYears(Data, Forc, FixedParams, ... 
    'singleYear', useSingleYear, 'singleSpectra', useSingleSizeSpectra);

% There are lots of forcing data particle trajectories -- far too many to
% include within a parameter optimisation procedure, so these need to be
% filtered.
% For each in-situ sample event (specified in Data.scalar), choose a set of
% particle trajectories that pass nearby.

% Choose trajectories from within maxDist km radius of sampling sites.
if ~exist('maxDist', 'var')
    maxDist = 25;
end
% Number of trajectories per sample event (duplicates may be required for 
% some events, depending on maxDist).
if ~exist('numTraj', 'var')
    numTraj = 10;
end

[Forc_, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, numTraj);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
Data = omitUnmatchedEvents(Data, eventTraj, Forc);
% [Data, eventTraj] = omitUnmatchedEvents(Data, eventTraj, Forc);

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
if ~exist('trajectoryPlot', 'var')
    trajectoryPlot = false;
end
if ~exist('dendrogramPlot', 'var')
    dendrogramPlot = false;
end

[Forc, Data.scalar] = particleOrigin(Forc, Data.scalar, ...
    'trajectoryPlot', trajectoryPlot, 'dendrogramPlot', dendrogramPlot); pause(0.25)

% Standardise the fitting data using linear mixed models to adjust for
% variability due to depth and sampling event.
if ~exist('plotScalarData', 'var')
    plotScalarData = false;
end
if ~exist('plotSizeData', 'var')
    plotSizeData = false;
end
if ~exist('plotAllData', 'var')
    plotAllData = false;
end

Data = standardiseFittingData(Data, ...
    'plotScalarData', plotScalarData, 'plotSizeData', plotSizeData, 'plotAllData', plotAllData);

% Include extra fields indexing sorted order of data -- convenience for
% plotting
Data.scalar = sortOrderData(Data.scalar);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Choose ODE solver
integratorChoices = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s'};
FixedParams.odeIntegrator = integratorChoices{2};

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ode23 seems to be the most robust. ode45 has occasionally produced NaNs.
% Maybe ode45 is less robust due to stiffness in the model equations. It
% seems that the equations are stiff for some parameter values, as the
% solver can slow down significantly... this is annoying, it would be
% useful to be able to dynamically switch between stiff and non-stiff
% solvers during the integrations (such switching functionality is
% available in the DEsolve package in R...)
% The stiff-solvers (ode15s, ode23s) are not working... don't know why...
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Set solver options.
% Integrator is called separately for each time step -- defined by forcing 
% data (daily) intervals.
odeMaxTime = FixedParams.dt_max; % max timestep (days)
odeInitTime = 0.5 * odeMaxTime; % initial integration timestep (solver will
                                % automatically reduce this if required)

FixedParams.odeOptions=odeset('AbsTol', 1e-6, 'RelTol', 1e-4,...
    'InitialStep', odeInitTime, 'MaxStep', odeMaxTime);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ~exist('displayAllOutputs', 'var')
    displayAllOutputs = false;
end
switch displayAllOutputs, case true
    displayForc = true;
    displayData = true;
    displayFixedParams = true;
    displayParams = true;
end

if ~exist('displayForc', 'var')
    displayForc = false;
end
switch displayForc, case true
    fprintf('\n\n')
    disp('Forc contains "particle trajectory" data from a physical model, which drives the size-structured NPZD model.')
    fprintf('\n\n')
    disp(['Forcing data dimensions: time = ' num2str(FixedParams.nt)])
    disp(['                         depth = ' num2str(FixedParams.nz)])
    disp(['                         trajectory = ' num2str(length(Forc.iTraj))])
    fprintf('\n\n')
    display(Forc)
end

if ~exist('displayData', 'var')
    displayData = false;
end
switch displayData, case true
    fprintf('\n\n')
    disp('Data contains the in-situ measurements used to optimise model parameters -- scalar and size-spectra data aggregated into separate structs.')
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

if ~exist('displayFixedParams', 'var')
    displayFixedParams = false;
end
switch displayFixedParams, case true
    fprintf('\n\n')
    disp('FixedParams contains all parameters with constant values, model dimensions, indexes etc -- parameters that cannot be numerically optimised.')
    display(FixedParams)
end

if ~exist('displayParams', 'var')
    displayParams = false;
end
switch displayParams, case true
    fprintf('\n\n')
    disp('Params contains model parameters that may be numerically optimised, and their bounding intervals.')
    display(Params)
end

end

