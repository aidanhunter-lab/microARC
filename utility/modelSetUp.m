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

% Chlorophyll samples from deep water cannot be replicated by the model
% unless modelled plankton sink -- light does not penetrate deeply enough 
% to permit growth and diffusion is too weak. Unless plankton sinking is
% modelled then the deep water chlorophyll samples should not be used when
% fitting model parameters.
ds = Data.scalar;
ds = rmfield(ds, {'nSamples','nEvents'});
ds = struct2table(ds);
ds(strcmp(ds.Variable, 'chl_a') & ds.Depth >= 100,:) = [];
ds = table2struct(ds, 'ToScalar', true);
ds.nSamples = length(ds.Value);
ds.nEvents = length(unique(ds.Event));
Data.scalar = ds;
clear ds

% Choose number of modelled size class intervals using function
% setSizeClassIntervals.m
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
if ~exist('nsizes', 'var')
   nsizes = 9;
end
% And cell diameter of smallest class
if ~exist('ESD1', 'var')
    ESD1 = 1.5;
end


% Initialise model parameters.
% Values can be modified in defaultParameters.m, which is called from
% within initialiseParameters.m.

% A file-path for saved parameters may be specified in setDirectories.m
% If a pre-existing parameter set is not loaded then set initialParamDir = [].
if ~exist('initialParamDir', 'var')
    initialParamDir = [];
end
[FixedParams, Params] = initialiseParameters(F, bioModel, ... 
    'ESDmin', ESDmin, 'ESDmax', ESDmax, 'ESD1', ESD1, 'nsizes', nsizes,  ...
    'parFile', initialParamDir);
nsizes = FixedParams.nPP_size;

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
% An extra few useful metrics are also calculated here.
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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% THIS CODE SECTION IS NOT ROBUST SO ANY CHANGES TO WHICH YEARS OF DATA ARE
% USED WILL NEED TO BE MADE WITH CARE, BEARING IN MIND THAT SUBSEQUENT CODE
% SECTIONS WILL BE AFFECTED, E.G., THE COST FUNCTION...

% Select which year(s) to model and filter out unused data.
% Function selectYears.m automatically chooses which year(s) to use based
% upon which years are most replete with data.

% Run model over a single year? If false, then multiple years of forcing 
% data MAY be used depending on fitting-data availability
if ~exist('useSingleYear', 'var'), useSingleYear = true; end
% In 2018 there are size spectra from 2 cruises --
% if useSingleSizeSpectra = true then use only the Polarstern samples
if ~exist('useSingleSizeSpectra', 'var'), useSingleSizeSpectra = true; end

[Data, Forc, FixedParams] = selectYears(Data, Forc, FixedParams, ... 
    'singleYear', useSingleYear, 'singleSpectra', useSingleSizeSpectra);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% There are lots of forcing data particle trajectories -- far too many to
% include within a parameter optimisation procedure, so these need to be
% filtered.
% For each in-situ sample event (specified in Data.scalar), choose a set of
% particle trajectories that pass nearby.

% Choose trajectories from within maxDist km radius of sampling sites.
if ~exist('maxDist', 'var')
    maxDist = 25;
end
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

% Number of trajectories per sample event (duplicates may be required for 
% some events, depending on maxDist).
if exist('numTraj', 'var')
    % numTrajDefault is used, in the first instance, to filter particle
    % trajectories while retaining sufficient (at least 10) trajectories
    % per sample event to determine the water mass origin of each sample.
    % Then, after water mass origin (Arctic/Atlantic/Arctic&Atlantic) is
    % determined, we re-filter to retain numTraj trajectories per sample.
    % This roundabout method is important because if numTraj=1 (which is
    % good for parameter optimisation efficiency) then samples at the
    % water mass boundaries are described as EITHER Arctic OR Atlantic,
    % which is dodgy because a sample may, in reality, be from Arctic water
    % while the nearest trajectory originated from the Atlantic, and we
    % cannot tell which samples are boundary cases from single trajecories.
    numTrajDefault = max(eval('numTraj'), 10);
else
    numTraj = 10;
    numTrajDefault = 10;
end

[Forc_, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, numTrajDefault, ...
    'displayWarnings', false);
% Observe any warnings produced by chooseTrajectories.m and use the table
% 'eventTraj' to help evaluate choice of maxDist and numTraj.

% Omit data from any sampling events not matched with a set of trajectories
Data_ = omitUnmatchedEvents(Data, eventTraj, Forc_);

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

[Forc_, Data_] = particleOrigin(Forc_, Data_, ...
    'trajectoryPlot', trajectoryPlot, 'dendrogramPlot', dendrogramPlot); pause(0.25)

% Now that water mass origin has been determined for each sample (using
% particle trajectories) re-filter to select numTraj particle trajectories
% per sample.
[Forc, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, numTraj);

% Omit data from any sampling events not matched with a set of trajectories
Data = omitUnmatchedEvents(Data, eventTraj, Forc);

% This code section is a bit messy because it was a quick-fix to sort out
% issues with ascribing water mass origins to sample events when using a
% single particle trajctory per event -- it could be improved...
Data.scalar.waterMass = Data_.scalar.waterMass;
Data.scalar.AtlanticOrigin = Data_.scalar.AtlanticOrigin;
Data.scalar.ArcticOrigin = Data_.scalar.ArcticOrigin;

Data.sizeFull.waterMass = Data_.sizeFull.waterMass;
Data.sizeFull.AtlanticOrigin = Data_.sizeFull.AtlanticOrigin;
Data.sizeFull.ArcticOrigin = Data_.sizeFull.ArcticOrigin;

Data.sizeFull.dataBinned.waterMass = Data_.sizeFull.dataBinned.waterMass;
Data.sizeFull.dataBinned.AtlanticOrigin = Data_.sizeFull.dataBinned.AtlanticOrigin;
Data.sizeFull.dataBinned.ArcticOrigin = Data_.sizeFull.dataBinned.ArcticOrigin;

Forc.waterMass = Forc_.waterMass(ismember(Forc_.iTraj, Forc.iTraj));

clear Forc_ Data_

% % Number of trajectories per sample event (duplicates may be required for 
% % some events, depending on maxDist).
% if ~exist('numTraj', 'var')
%     numTraj = 10;
% end
% 
% [Forc_, eventTraj] = chooseTrajectories(Forc, Data.scalar, maxDist, numTraj);
% 
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % The choice of maxDist is important. If maxDist is very large then the 
% % trajectories ascribed to each sampling event will correspond poorly to
% % the history of the water at the events, so small values are better.
% % However, if maxDist is small then this can limit the number of
% % trajectories available to use for each sampling event, and may entirely
% % exclude some sampling events from the analyses if no trajectories are
% % within the maxDist radius. We want to set maxDist small, but not too
% % small...
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% % Observe any warnings produced by chooseTrajectories.m and use the table
% % 'eventTraj' to help evaluate choice of maxDist and numTraj: if happy then
% % set Forc = Forc_, otherwise reselect maxDist and numTraj and repeat
% % chooseTrajectories.m
% Forc = Forc_; clear Forc_
% 
% % Omit data from any sampling events not matched with a set of trajectories
% Data = omitUnmatchedEvents(Data, eventTraj, Forc);
% % [Data, eventTraj] = omitUnmatchedEvents(Data, eventTraj, Forc);
% 
% % Group sampling events by origin of particles: Arctic or Atlantic.
% % Each individual trajectory is either of Atlantic or Arctic origin
% % (stored in Forc.waterMass).
% % Each sampling event is associated with a selection of trajectories, so
% % may be considered as Atlantic, Arctic, or a mixture (stored in 
% % Data.scalar.waterMass).
% if ~exist('trajectoryPlot', 'var')
%     trajectoryPlot = false;
% end
% if ~exist('dendrogramPlot', 'var')
%     dendrogramPlot = false;
% end
% 
% [Forc, Data] = particleOrigin(Forc, Data, ...
%     'trajectoryPlot', trajectoryPlot, 'dendrogramPlot', dendrogramPlot); pause(0.25)


% Group size data by water origin -- find average spectra using
% measurements from events whose trajectories all orginate from either the
% Arctic or the Atlantic
% avFun = @mean; % arithmetic or geometric mean spectra over sampe events and depths
avFun = @geomean; % choice of average type is important! geometric mean more appropriate for these data?
Data = sizeDataOrigin(Data, 'avFun', avFun);

% For each trajectory, find the time of the latest sampling event.
% Integrations along trajectories can then be stopped at these events to
% reduce model run-times during parameter optimisation.
Forc = latestSampleTime(Forc, Data);

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

if isfield(FixedParams,'NclineDepth')
    NclineDepth = FixedParams.NclineDepth;
else
    NclineDepth = [];
end

Data = standardiseFittingData(Data, ...
    'plotScalarData', plotScalarData, 'plotSizeData', plotSizeData, ...
    'plotAllData', plotAllData, 'NclineDepth', NclineDepth);

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
% Integration tolerences
if ~exist('RelTol', 'var')
    RelTol = 1e-2; % RelTol = 1e-2 => results accurate to 1%
end
if ~exist('AbsTol', 'var')
    AbsTol = 1e-3;
end

FixedParams.odeOptions=odeset('AbsTol', AbsTol, 'RelTol', RelTol,...
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

