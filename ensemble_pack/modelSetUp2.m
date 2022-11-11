% this is a manipulated version of modelSetUp.m
% sets up model for trajectories that cross sampling region in summer or autumn

% CHANGES
    % varargin seasonConfig specifies whether model is run for 'summer' 
    %   (default) or 'autumn'
    % load new prepareFittingData2 (loads all data, summer and autumn),
    % then subsetting happens


function [Forc, FixedParams, Params, Data] = modelSetUp2(Directories, varargin)

% modelSetup.m is a wrapper for numerous functions used to load and
% organise all necessary (forcing and fitting) data, and to initialise
% model parameters.

% Model input is stored in 4 structs: Forc, FixedParams, Params, and Data,
% returned by modelSetup.m
% Forc        = forcing data from physical model particle trajectories.
% FixedParams = constant model parameters (not numerically tuned).
% Params      = model parameters that MAY be tuned.
% Data        = observed data included in cost function to tune parameters.

% Plots may be generated as extra output by setting as true any of the
% following optional arguments coded in varargin as name-value pairs...
% 'plotCellConcSpectra'
% 'plotBioVolSpectra'
% 'plotSizeClassIntervals'
% 'trajectoryPlot'
% 'dendrogramPlot'
% 'plotScalarData'
% 'plotSizeData'
% 'plotAllData'

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

if ~exist('seasonConfig', 'var')
    seasonConfig = 'summer';
    warning("No seasonConfig selected. Using 'summer' as default.")
end


if ~any([strcmp(seasonConfig, 'summer') strcmp(seasonConfig, 'autumn')])
    error('seasonConfig must be either "summer" or "autumn".' )
end



% Load and extract relevant forcing data from physical model 'forcModel'.
switch seasonConfig
    case 'summer'
        % Specify years of available forcing data
        if ~exist('years', 'var')
%             years = 2017:2018;  % das gibt später probleme, denn
%             watermass origin kann immer nur für ein jahr gemacht werden.
%             Per default nehmen wir also 2018, aber 2017 (sommer) wäre
%             auch möglich
            years = 2018; 
            warning("'years' were not specified, using default setting (2018)")
        end

        % Load and extract relevant forcing data from physical model 'forcModel'.
        F = loadForcing(Directories, years, useTraj);
        % ~~~
        
    case 'autumn'
        % load Forcing for autumn trajectories, only available for 2018
        
        % Specify years of available forcing data
        if ~exist('years', 'var')
            years = 2018;
        else
            % warning if year is specified and not 2018, then break break
            if any(years ~= 2018)
                error("Forcing for autumn trajectories is only available for 2018. No Forcing loaded.")
                % return 
            end
        end
                   
        % assign new path to autumn forc file 
        Directories.forcDummy = 'particles_MODEL_DOMAIN_EXPERIMENT-YEARb_YEAR_t*iSub03.mat'; % das b in YEARb ist entscheidend
        F = loadForcing(Directories, years, useTraj);
end

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
Data = prepareFittingData2(obsDir, ...
    'plotCellConcSpectra', plotCellConcSpectra, ...
    'plotBioVolSpectra', plotBioVolSpectra);
% % Data now contains data from all 4 cruises. 


% subset data to only contain relevant season AND YEAR(S)
switch seasonConfig
    case 'summer'
        
        % subset scalar data
        fnames = fields(Data.scalar); 
        ind_season = ~contains(Data.scalar.Label, 'MSM77');
        ind_year = Data.scalar.Year ==  years; 
        ind = ind_season & ind_year;
        for i = 1:length(fnames)-2
            
            Data_.scalar.(fnames{i}) = Data.scalar.(fnames{i})(ind); 
        end
        Data_.scalar.nSamples = length(Data_.scalar.Label); 
        Data_.scalar.nEvents = length(unique(Data_.scalar.Event)); 
        clear fnames ind
        
        % subset size data
        fnames = fields(Data.size);
        ind = contains(Data.size.scenario, {'S1', 'S2'}); % mean of 2016, 2017 and 2018 specra
        
        for i = 1:length(fnames)
            Data_.size.(fnames{i}) = Data.size.(fnames{i})(ind); 
        end
        clear fnames ind
        
        % subset sizeFull data
        fnames = fields(Data.sizeFull);
        ind_season = ~contains(Data.sizeFull.Label, 'MSM77');
        ind_year = Data.sizeFull.Year == years; 
        ind = ind_season & ind_year; 
        
        for i = 1:length(fnames)
            Data_.sizeFull.(fnames{i}) = Data.sizeFull.(fnames{i})(ind);
        end
        clear fnames ind
        
        
        
    case 'autumn'  % there is only 2018 anyway
        
        % subset scalar data
        fnames = fields(Data.scalar); 
        ind = contains(Data.scalar.Label, 'MSM77');
        for i = 1:length(fnames)-2
            
            Data_.scalar.(fnames{i}) = Data.scalar.(fnames{i})(ind); 
        end
        Data_.scalar.nSamples = length(Data_.scalar.Label); 
        Data_.scalar.nEvents = length(unique(Data_.scalar.Event)); 
        clear fnames ind
        
        % subset size data
        fnames = fields(Data.size);
        ind = contains(Data.size.scenario, {'S3', 'S4'}); 
        
        for i = 1:length(fnames)
            Data_.size.(fnames{i}) = Data.size.(fnames{i})(ind); 
        end
        clear fnames ind
        
        % subset sizeFull data
        fnames = fields(Data.sizeFull);
        ind = contains(Data.sizeFull.Label, 'MSM77');
        
        for i = 1:length(fnames)
            Data_.sizeFull.(fnames{i}) = Data.sizeFull.(fnames{i})(ind);
        end
        clear fnames ind
    
end

Data = Data_ ; 

clear Data_
% relabel the event numbers, OR ELSE chooseTrajectories.m will fail.
    eventLabs = table((1:Data.scalar.nEvents)', unique(Data.scalar.Event));
    eventLabs.Properties.VariableNames = {'newEvent', 'Event'};
    tmp = table(Data.scalar.Event); tmp.Properties.VariableNames = {'Event'};
    tmp = join(tmp, eventLabs);
    Data.scalar.Event = tmp.newEvent;

    tmp = table(Data.sizeFull.Event); tmp.Properties.VariableNames = {'Event'};
    tmp = join(tmp, eventLabs);
    Data.sizeFull.Event = tmp.newEvent;

%     tmp = table(Data.sizeFull.dataBinned.Event); tmp.Properties.VariableNames = {'Event'};
%     tmp = join(tmp, eventLabs);
%     Data.sizeFull.dataBinned.Event = tmp.newEvent;




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

% Choose data types to use in cost function
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
% clear F

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% The 4 main structs (Forc, FixedParams, Params, Data) should now be stored
% in the workspace.
% The rest of this set-up code section selects which data to use, 
% standardises data, filters out unrequired data, and derives a few useful
% metrics.

% Remove fitting-data samples from below the maximum modelled depth.
% Optional arguments may be included to omit samples from below selected depths.

% Chlorophyll samples from deep water cannot be replicated by the model
% unless modelled plankton sink -- light does not penetrate deeply enough 
% to permit growth and diffusion is too weak. Unless plankton sinking is
% modelled then the deep water chlorophyll samples should probabaly not be 
% used when fitting model parameters.
if ~exist('chlSampleDepthLimit', 'var')
    chlSampleDepthLimit = 100; % By default exclude chl samples from depths >= 100
end

Data = omitDeepSamples(Data, FixedParams, ...
    'chlSampleDepthLimit', chlSampleDepthLimit);
Data.sizeFull.nSamples = length(Data.sizeFull.Label); % wird später gebraucht, passiert sonst eig in selectYears


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % % THIS CODE SECTION IS NOT ROBUST SO ANY CHANGES TO WHICH YEARS OF DATA ARE
    % % USED WILL NEED TO BE MADE WITH CARE, BEARING IN MIND THAT SUBSEQUENT CODE
    % % SECTIONS WILL BE AFFECTED, E.G., THE COST FUNCTION...
    % 
    % % Select which year(s) to model and filter out unused data.
    % % Function selectYears.m automatically chooses which year(s) to use based
    % % upon which years are most replete with data.
    % 
%     % Run model over a single year? If false, then multiple years of forcing 
%     % data MAY be used depending on fitting-data availability
%     if ~exist('useSingleYear', 'var'), useSingleYear = true; end
%     % In 2018 there are size spectra from 2 cruises --
%     % if useSingleSizeSpectra = true then use only the Polarstern samples
%     if ~exist('useSingleSizeSpectra', 'var'), useSingleSizeSpectra = true; end         %% diese varargin muss true sein! damit MSM77 mit drin bleibt, wenn es so sein soll. ODER?? Data ist ja eh gefiltert, und sollte nur PS114 ODER MSM77 haben...  
%     
%     [Data, Forc, FixedParams] = selectYears(Data, Forc, FixedParams, ... 
%         'singleYear', useSingleYear, 'singleSpectra', useSingleSizeSpectra);
%     
    
%         % Run model over a single year? If false, then multiple years of forcing 
%     % data MAY be used depending on fitting-data availability
%     if ~exist('useSingleYear', 'var'), useSingleYear = true; end
%     % In 2018 there are size spectra from 2 cruises --
%     % if useSingleSizeSpectra = true then use only the Polarstern samples
%     if ~exist('useSingleSizeSpectra', 'var'), useSingleSizeSpectra = true; end         %% diese varargin muss true sein! damit MSM77 mit drin bleibt, wenn es so sein soll. ODER?? Data ist ja eh gefiltert, und sollte nur PS114 ODER MSM77 haben...  
%     
%     [Data_, Forc_, FixedParams_] = selectYears(Data, Forc, FixedParams, ... 
%         'singleYear', useSingleYear, 'singleSpectra', useSingleSizeSpectra); %%% SOLL ein test sein, funktioniert aber nicht. nichts funktioniert hier. 
 % Das dunktioniert nicht, weil obwhol year = 2018 gestetzt ist, das jahr
 % mit den meisten Daten gewählt wird. für spektren ist das 2017... daher
 % stimmen die events nicht mehr überein und es gibt 1 error....
    
    % % THE COST FUNCTIONS WONT WORK ANYMORE, ADJUST matchModOutput2Data!!
    
    % % AUCH: der part hierdrunter, watermass origin zuornung, muss
    % repariert werden. 


    
    % new approach: filter Data for year; so that trajectory slection
    % further downstream still works...
    
    
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% If forcing data have not already been grouped by particle origin (Arctic
% or Atlantic) then assign each particle to a group. This is time-consuming
% when using lots of particles, so the function is run once and the output
% saved.

for ii = 1:length(Forc.years)
    year = Forc.years(ii);
    
    switch seasonConfig
        case 'summer'
            fileName = ['particleGroupings_' num2str(year) '.mat'];
        case 'autumn' 
            fileName = ['particleGroupings_' num2str(year) 'b.mat'];
    end
    loadGroups = exist(fileName, 'file');
    if loadGroups
        loadedGroups = load(fullfile(forcDir,fileName));
        particleGrouping = loadedGroups.particleGrouping;
    end
    match = length(Forc.iTraj) == height(particleGrouping); % do loaded values match forcing data?
    
    if loadGroups && match
        Forc.waterMass = particleGrouping.waterMass';
    else
        warning('The particleOrigin.m function has been called to group all forcing data particles as either Arctic or Atlantic origin. This time-consuming function should only need called once per set of forcing data as the outputs are saved and re-used. This warning should only appear the 1st time modelSetUp.m is called for a new set of forcing data, otherwise there may be some problem.')
        
        if ~exist('trajectoryPlot', 'var')
            trajectoryPlot = false;
        end
        if ~exist('dendrogramPlot', 'var')
            dendrogramPlot = false;
        end

        Forc = particleOrigin(Forc, 'trajectoryPlot', trajectoryPlot, ... 
            'dendrogramPlot', dendrogramPlot, 'progressBar', true); pause(0.25)
        
        % Extract and save water mass groupings so they may simply be
        % loaded next time modelSetUp.m is called.
        
        particleGrouping = table(repmat(year, [Forc.nTraj, 1]), (1:Forc.nTraj)', Forc.waterMass(:), ...
            'VariableNames', {'year','traj','waterMass'});
        
        save(fullfile(forcDir, fileName), 'particleGrouping')
        
    end

end

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

[~, Data_] = particleOrigin(Forc_, 'Data', Data_, ...
    'trajectoryPlot', trajectoryPlot, 'dendrogramPlot', dendrogramPlot); pause(0.25)
% [Forc_, Data_] = particleOrigin(Forc_, 'Data', Data_, ...
%     'trajectoryPlot', trajectoryPlot, 'dendrogramPlot', dendrogramPlot); pause(0.25)

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

% Forc.waterMass = Forc_.waterMass(ismember(Forc_.iTraj, Forc.iTraj));

clear Forc_ Data_


% Group size data by water origin -- find average spectra using
% measurements from events whose trajectories all orginate from either the
% Arctic or the Atlantic
if ~exist('avFun_sizeSpectra', 'var')
    avFun_sizeSpectra = @geomean; % by default find averages of multiple size spectra using geometric means
end

Data = sizeDataOrigin(Data, 'avFun', avFun_sizeSpectra);

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

if isfield(FixedParams,'NclineDepth') && strcmp(seasonConfig, 'summer')    % add summmer condition here, so that NclineDepth = [] when 'autumn' and Data standardisation can work
    NclineDepth = FixedParams.NclineDepth;
else
    NclineDepth = [];
end

% data standardisatin only works as is for the 'summer' case.
% for the standardiseFittingData function to work without further
% manipulation for the 'autumn' case, we first create dummy entries for
% the missing variables (N, PON) that will be removed again after the
% process.

if strcmp(seasonConfig, 'autumn')
    
    % data types needed for standardiseFittingData function:
    obsRequired = Data.scalar.obsInCostFunction; 
    % data types that are available: 
    obsAvailable = unique(Data.scalar.Variable);

    % which obs are missing?
    [shared, index] = intersect(obsRequired, obsAvailable);
    obsMissing = obsRequired;
    obsMissing(index) = [];
    clear shared index obsAvailable obsRequired
    
    
    temp = Data.scalar; 
    fnames = fieldnames(temp);
    nSamples_orig = Data.scalar.nSamples;
    % attach dummies for .Variable
    temp.Variable = vertcat(temp.Variable, repelem(obsMissing', 2)); 
    %attach dummies for all other nSample-long fields, so table can be built in standardiseFittingData
    for i = 1:length(fnames)
        if size(temp.(fnames{i}),1) == temp.nSamples 
            temp.(fnames{i}) = vertcat(temp.(fnames{i}), repmat(temp.(fnames{i})(end-1:end),[2 1]));    
        end
    end
    % and change nSamples
    temp.nSamples = temp.nSamples + 2*length(obsMissing);
    % reaasign temp to Data.scalar
    Data.scalar = temp;
    clear temp obsMissing
end


Data = standardiseFittingData(Data, ...
    'plotScalarData', plotScalarData, 'plotSizeData', plotSizeData, ...
    'plotAllData', plotAllData, 'NclineDepth', NclineDepth);

% if 'autumn', remove the dummy entries again. 
if strcmp(seasonConfig, 'autumn')
    Data.scalar.nSamples = nSamples_orig;
    for i = 1:length(fnames)
        if size(Data.scalar.(fnames{i}),1) >= Data.scalar.nSamples 
            Data.scalar.(fnames{i}) = Data.scalar.(fnames{i})(1:Data.scalar.nSamples);    
        end
    end 
    clear nSamples_orig fnames
end


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

