function [Data, eventTraj] = omitUnmatchedEvents(Data, eventTraj, Forc)
% Omit data from events that have not been ascribed to trajectories

uevents = unique(eventTraj.event);
nevents = length(uevents);

% Filter data
if Data.scalar.nEvents > nevents
    Data.scalar.nEvents = nevents;
    keep = ismember(Data.scalar.Event, uevents);
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

% Match sampling events to trajectories
Data.scalar.EventTraj = zeros(nevents, Forc.nTraj);
for i = 1:nevents
    ti = eventTraj.trajIndex(eventTraj.event == i);
    for j = 1:length(ti)
        Data.scalar.EventTraj(i,ti(j)) = Data.scalar.EventTraj(i,ti(j)) + 1;
    end
end
for i = 1:nevents
    tn = Data.scalar.EventTraj(i,:); % trajectories -- a few may be duplicates
    ntraj = sum(tn);
    if i == 1
        Data.scalar.evTraj = nan(ntraj, nevents);
    end
    ti = find(tn > 0);
    % a few events may have duplicate trajectories...
    indrep = tn(ti);
    ti = arrayfun(@(x) repmat(ti(x), indrep(x), 1), 1:numel(ti), 'uni', 0);
    Data.scalar.evTraj(:,i) = vertcat(ti{:});
end

% Rearrange the trajectories to minimise duplicates across the samples... 
% this is beneficial for the cost function
nsamples = size(Data.scalar.evTraj, 1);
for i = 2:size(Data.scalar.evTraj, 2) % this loop is not robust code, but should be OK if sample number is not very large
    xp = Data.scalar.evTraj(:,1:i-1);
    xi = Data.scalar.evTraj(:,i);    
    duplicateTrajectory = repmat(xi, [1, i-1]) == xp;
    duplicateTrajectory = any(duplicateTrajectory(:));    
    j = 1;
    while duplicateTrajectory
        xi = circshift(xi, 1);
        duplicateTrajectory = repmat(xi, [1, i-1]) == xp;
        duplicateTrajectory = any(duplicateTrajectory(:));
        j = j+1;
        if j == nsamples
            xi = Data.scalar.evTraj(:,i);
            break;
        end
    end
    Data.scalar.evTraj(:,i) = xi;
    eventTraj.trajIndex(eventTraj.event == i) = xi;
end



% Size data
discard = ~ismember(Data.sizeFull.Label, Data.scalar.Label);
if any(discard)
    fields = fieldnames(Data.sizeFull);
    for i = 1:length(fields)
        if size(Data.sizeFull.(fields{i}), 1) == length(discard)
            Data.sizeFull.(fields{i}) = Data.sizeFull.(fields{i})(~discard);
        end
    end
    Data.sizeFull.nSamples = length(Data.sizeFull.Label);
end

tt = unique(table(Data.scalar.Label, Data.scalar.Event));
tt.Properties.VariableNames = {'Label','Event'};

sizeData = Data.sizeFull;
sizeData = rmfield(sizeData, {'dataBinned','obsInCostFunction','nSamples'});
sizeData = struct2table(sizeData);

sizeData.Event = [];
sizeData = innerjoin(sizeData, tt);

sizeData = movevars(sizeData, 'Event', 'After', 'Sample');

sizeData = table2struct(sizeData, 'ToScalar', true);

sizeData.dataBinned = Data.sizeFull.dataBinned;
sizeData.obsInCostFunction = Data.sizeFull.obsInCostFunction;
sizeData.nSamples = Data.sizeFull.nSamples;

Data.sizeFull = sizeData;


sizeData = struct2table(Data.sizeFull.dataBinned);

sizeData.Event = [];
sizeData = innerjoin(sizeData, tt);
sizeData = movevars(sizeData, 'Event', 'After', 'Yearday');
sizeData = table2struct(sizeData, 'ToScalar', true);

Data.sizeFull.dataBinned = sizeData;

