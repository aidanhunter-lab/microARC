function [FixedParams, Params, Forc, Data] = optimisationOptions(FixedParams, Params, Forc, Data, varargin)
% Choose tuning parameters and cost function and numerical tuning algorithm
% and any other options related to optimisation can be included here...

extractVarargin(varargin)

%% Select parameters to optimise

% Choose from the lists: Params.scalars & Params.sizeDependent.
parnames = {'wPOM1', 'wp_a', 'wp_b', 'rDON', 'rPON', ...
    'aP', 'm_b', 'Gmax_a', 'Gmax_b', 'k_G', 'pmax_a', 'pmax_b', ... 
    'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', 'Qmax_delQ_b', ... 
    'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b'};

% Check that all chosen parameters exist -- error if not
if ~all(ismember(parnames, Params.scalarParams) | ... 
        ismember(parnames, Params.vectorParams))
    error('Error in "optimisationOptions.m". Invalid choice for "parnames": all tuning parameter names must appear in "Params.scalarParams" or "Params.vectorParams".')
end

% Extract tuning parameter initial values and bounds
npars = length(parnames);
par0 = cell2mat(cellfun(@(x) Params.(x), parnames, 'UniformOutput', false));
bounds = cellfun(@(x) Params.bounds.(x), parnames, 'UniformOutput', false);
lb = cellfun(@(x) x(1), bounds); % lower bounds
ub = cellfun(@(x) x(2), bounds); % upper bounds
% Store tuning parameter names and their bounding values
FixedParams.tunePars = parnames;
FixedParams.tunePars_lb = lb;
FixedParams.tunePars_ub = ub;
% Assign to workspace
assignin('caller', 'npars', npars)
assignin('caller', 'boundsLower', lb)
assignin('caller', 'boundsUpper', ub)

%% Select cost function.
% There's a few options for the cost function. Not yet sure which is the
% best... Hellinger2_groupWaterOrigin is most defensible as it makes fewest
% assumptions
[~, ~, costFunctionChoices] = costFunction();
% costFunctionChoices = { ...
%     'LSS', ...
%     'RMS', ...
%     'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal', ...
%     'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet', ...
%     'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet', ...
%     'N_LN-Dir_groupWaterOrigin', ...
%     'Hellinger_groupWaterOrigin', ...
%     'Hellinger2_groupWaterOrigin', ...
%     'Hellinger_MVN_groupWaterOrigin'
%     };

if ~exist('costFunctionType', 'var')
    % costFunctionChoices should be given shorter names...
%     costFunctionType = costFunctionChoices{4}; % select cost function
    costFunctionType = costFunctionChoices{6}; % select cost function
end
FixedParams.costFunction = costFunctionType;
assignin('caller', 'costFunctionLabel', costFunctionType)

%% Select optimising algorithm (so far only ga is available)
optimiserChoices = {'ga','muga'};
if ~exist('optimiser', 'var')
    optimiser = optimiserChoices{1};
end
FixedParams.optimiser = optimiser;
optimise = str2func(optimiser);
assignin('caller', 'optimise', optimise) % assign optimising algorithm to the workspace


%% Parameters of optimising algorithm

% Parameters are optimised using a numerical population-based algorithm
if ~exist('popSize', 'var')
    popSize = 100; % number of parameter sets in algorithm population    
end
if ~exist('niter', 'var')
    niter = 10; % algorithm iterations
end
assignin('caller', 'popSize', popSize)
assignin('caller', 'niter', niter)

% Optimising algorithm options
switch optimiser
    case 'ga'
        optimiserOptions = optimoptions('ga', ...
            'PopulationSize', popSize, ...
            'Generations', niter, ...
            'InitialPopulationRange', [lb;ub], ...
            'InitialPopulationMatrix', par0, ...
            'Display', 'iter', ...
            'OutputFcn', @gaStoreHistory, ... % gaStoreHistory.m function dynamically stores output, which is available in workspace even if algorithm is halted prematurely
            'PlotFcn', @gaplotbestf);
end
assignin('caller', 'optimiserOptions', optimiserOptions)

% Halt integrations at each trajectory's final sampling events to save time
Forc.integrateFullTrajectory = false; % Integrating full trajectories not required for optimisation


%% Filter forcing- and fitting-data

% The default set-up uses all trajectories from Arctic and Atlantic waters
% and fits to the size data that's aggregated over all sample events
% (Data.size).
default = ~exist('fitTrajectories', 'var') || isempty(fitTrajectories); % for default set-up do not specify, or set as empty, fitTrajectories

% We can fit to size data that are grouped as Arctic or Atlantic samples.
% This allows more precision in matching in-situ samples to particle trajectories.
% Setting fitTrajectories = 'Atlantic', 'Arctic' or 'all' will use data
% derived from samples in that region and discard unneccessary data.

if ~default
    waterMasses = unique(Forc.waterMass);
    if ~ismember(fitTrajectories, [waterMasses, 'all'])
        error('Optional name-value pair argument "fitTrajectories" must be "Atlantic", "Arctic" or "all".')
    end
    
    % Filter forcing data - discard trajectories not used to fit model
    switch fitTrajectories
        case waterMasses
            keepTraj = strcmp(fitTrajectories, Forc.waterMass); % trajectories to retain
        case 'all'
            keepTraj = true(size(Forc.waterMass));
    end
    nkeep = sum(keepTraj); % updated number of trajectories after filtering
    ntraj = length(Forc.iTraj);
    fields = fieldnames(Forc);
    for i = 1:length(fields)
        x = Forc.(fields{i});
        nd = ndims(x);
        filter = size(x, nd) == ntraj; % trajectories stored in last dimension
        if filter
            x = shiftdim(x, nd-1);
            xs = size(x);
            x = x(keepTraj,:);
            x = reshape(x, [nkeep, xs(2:nd)]);
            x = shiftdim(x, 1);
            Forc.(fields{i}) = x;
        end
    end
    Forc.nTraj = length(Forc.iTraj);
    
    % Create table to relabel trajectories if some have been omitted
    oldTraj = (1:ntraj)';
    newTraj = keepTraj(:);
    relabelTraj = table(oldTraj, newTraj);
    relabelTraj = relabelTraj(relabelTraj.newTraj,:);
    relabelTraj.newTraj = (1:height(relabelTraj))';
        
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Filter fitting data - discard data from sample events associated with
    % the discarded trajectories.
    
    % Scalar data
    switch fitTrajectories
        case waterMasses
            keepEvents = strcmp(fitTrajectories, Data.scalar.waterMass); % index retained data by sample event
        case 'all'
            keepEvents = true(size(Data.scalar.waterMass));
    end
    Data.scalar.EventTraj = Data.scalar.EventTraj(keepEvents,keepTraj);
    Data.scalar.evTraj = Data.scalar.evTraj(:,keepEvents);
    Data.scalar.evTraj = arrayfun(@(x) ... 
        relabelTraj.newTraj(relabelTraj.oldTraj == x) , Data.scalar.evTraj); % relabel trajectories to match the number stored in Forc
    Data.scalar.waterMass = Data.scalar.waterMass(keepEvents);
    Data.scalar.AtlanticOrigin = Data.scalar.AtlanticOrigin(keepEvents);
    Data.scalar.ArcticOrigin = Data.scalar.ArcticOrigin(keepEvents);
    Data.scalar.nEvents = sum(keepEvents);
    keepData = ismember(Data.scalar.Event, find(keepEvents)); % index retained data
    nkeep = sum(keepData);
    fields = fieldnames(Data.scalar);
    nSamples = Data.scalar.nSamples;
    for i = 1:length(fields)
        x = Data.scalar.(fields{i});
        xs = size(x);
        nd = ndims(x);
        filter = xs(1) == nSamples;
        if filter
            x = x(keepData,:);
            x = reshape(x, [nkeep, xs(2:nd)]);
            Data.scalar.(fields{i}) = x;
        end
    end
    Data.scalar.nSamples = length(Data.scalar.t);
    Data.scalar = sortOrderData(Data.scalar);
    
    % Size data
    sizeData = Data.sizeFull;
    ev = unique(sizeData.Event);
    switch fitTrajectories
        case waterMasses
            keepEvents = strcmp(fitTrajectories, sizeData.waterMass); % index retained events
        case 'all'
            % Keep all samples asociated with either Arctic OR Atlantic
            % water -- discard data associated with a mix of Arctic and 
            % Atlantic trajectories. Hence we use separate size data for
            % trajectories originating from each region, and omit samples
            % near the boundaries because their origin is uncertain.
            keepEvents = ismember(sizeData.waterMass, waterMasses);
    end
    sizeData.EventTraj = sizeData.EventTraj(keepEvents,keepTraj);
    sizeData.evTraj = sizeData.evTraj(:,keepEvents);
    sizeData.evTraj = arrayfun(@(x) ... 
        relabelTraj.newTraj(relabelTraj.oldTraj == x) , sizeData.evTraj); % relabel trajectories to match the number stored in Forc
    sizeData.waterMass = sizeData.waterMass(keepEvents);
    sizeData.AtlanticOrigin = sizeData.AtlanticOrigin(keepEvents);
    sizeData.ArcticOrigin = sizeData.ArcticOrigin(keepEvents);
    keepData = ismember(sizeData.Event, ev(keepEvents)); % index retained data
    nSamples = sizeData.nSamples;
    fields = fieldnames(sizeData);
    for i = 1:length(fields)
        x = sizeData.(fields{i});
        xs = size(x);
        filter = xs(1) == nSamples;
        if filter
            x = x(keepData);
            sizeData.(fields{i}) = x;
        end
    end
    sizeData.nSamples = length(sizeData.t);
    
    % Binned data
    dataBinned = sizeData.dataBinned.groupedByOrigin;
    switch fitTrajectories
        case waterMasses
            keepData = strcmp(fitTrajectories, dataBinned.waterMass); % index retained data
        case 'all'
            keepData = ismember(dataBinned.waterMass, waterMasses);
    end
    dataBinned = structfun(@(x) x(keepData), dataBinned, ...
        'UniformOutput', false); % filter
    sizeData.dataBinned = dataBinned;
    
    % Simplify the size data struct
    Data = rmfield(Data, {'size','sizeFull'});
    Data.size = sizeData;
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % Renumber sample events -- starting at 1
    eventsOld = Data.scalar.Event(:);
    et = table(eventsOld);
    eventsOld = unique(eventsOld);
    eventsNew = (1:length(eventsOld))';
    et_ = table(eventsOld, eventsNew);
    et = innerjoin(et, et_);
    Data.scalar.Event = et.eventsNew;
    
    eventsOld = Data.size.Event(:);
    et = table(eventsOld);
    et = innerjoin(et, et_);
    Data.size.Event = et.eventsNew;
    
else
    % Default - fitTrajectories either empty or not specified => use all
    % particles trajectories, and the size data that's aggregated over all
    % samples.
    
    % Remove unneccessary data
    Data = rmfield(Data, 'sizeFull');
    
end


