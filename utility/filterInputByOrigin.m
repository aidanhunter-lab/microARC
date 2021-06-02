function [Forc, Data] = filterInputByOrigin(Forc, Data, varargin)

% Model input structs containing forcing and fitting data are filtered
% according to the origin of the forcing data particle trajectories.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% An optional input argument, 'fitTrajectories', determines how data are
% filtered...
%
% fitTrajectories = 'Arctic' or 'Atlantic': forcing data filtered to only
%                                           include trajectories originating
%                                           from selected region; fitting data 
%                                           also filtered to only include
%                                           samples associated with the
%                                           retained trajectories.
% fitTrajectories = 'all': data are not filtered; separate size data used 
%                          for Arctic and Atlantic trajectories... fitting
%                          to Arctic and Atlantic size data simultaneously
%                          is problematic => avoid this option when
%                          optimising model parameters.
% fitTrajectories = [] or not specified: data are not filtered; single size
%                                        data vectors used for Arctic and
%                                        Atlantic trajectories -- size data
%                                        used for fitting are fully aggregated
%                                        over all sampling events.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extractVarargin(varargin)

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
%     Data.scalar = sortOrderData(Data.scalar);

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
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % If data, not needed for optimisation, has been discarded then
    % re-standardise the remaining data
    Data = standardiseFittingData(Data);
    Data.scalar = sortOrderData(Data.scalar);

    
else
    % Default - fitTrajectories either empty or not specified => use all
    % particles trajectories, and the size data that's aggregated over all
    % samples.
    
    % Remove unneccessary data
    Data = rmfield(Data, 'sizeFull');
    
end






