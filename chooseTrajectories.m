function [Forc, events] = chooseTrajectories(Forc, dat, maxDist, maxTraj)

% Select trajectories for each sampling event based on distance between
% particles and sample.
% Filter out any unused forcing data.
% Return a table linking event numbers to particle indexes.

nt = size(Forc.t,1);
nTraj = length(Forc.iTraj);

% particle times as year and yearday
time = yearday(Forc.t);
year = datevec(Forc.t);
year = reshape(year(:,1), [nt nTraj]);

% Find distance between particles and event at event time
nEvent = length(unique(dat.event));
Dist = nan(nTraj, nEvent);
Rank = nan(size(Dist));
Keep = nan(size(Dist));

for i = 1:nEvent
    ie = dat.event == i;
    sYear = dat.year(ie);
    sLat = dat.lat(ie);
    sLon = dat.lon(ie);
    sTime = dat.yearday(ie);
    samplePos_rads = deg2rad([sLon sLat]);

    keep2d = sTime == time & sYear == year;    
    keep = any(keep2d);
    
    particlePos_rads = deg2rad([Forc.x(keep2d) Forc.y(keep2d)]);
    Dist(keep,i) = distance_on_earth(samplePos_rads, particlePos_rads) / 1000;
    
    % Filter by maximum distance from event
    nearby = Dist(:,i) <= maxDist;
    keep = keep(:) & nearby(:);
    Dist(~keep,i) = nan;
    
    % Rank particles by distance
    [~,I] = sort(Dist(keep,i));
    [~,rankDist] = sort(I);    
    Rank(keep,i) = rankDist;
    
    % Limit number of trajectories to maxTraj per event
    nbest = Rank(:,i) <= maxTraj;
    keep = keep(:) & nbest(:);
    Dist(~keep,i) = nan;
    Rank(~keep,i) = nan;
    
    Keep(:,i) = keep;
    
end


iTraj = any(Keep,2); % index selected trajectories
nTraj_new = sum(iTraj);
Forc.nTraj = nTraj_new;
fields = fieldnames(Forc);
nfields = length(fields);
Forc = orderfields(Forc, [1 nfields 2:nfields-1]);

trajIndex = find(iTraj);
trajIndex_new = 1:nTraj_new;
trajIndex = trajIndex(:); trajIndex_new = trajIndex_new(:);
indexChange = table(trajIndex, trajIndex_new); % filtering trajectories changes the indexing

% store event numbers and associated trajectories in a table
nlinks_e = sum(Keep);
nlinks = sum(nlinks_e);
events = table(nan(nlinks,1),nan(nlinks,1));
events.Properties.VariableNames = {'event', 'trajIndex'};
j=1; k=0;
for i = 1:nEvent    
    ind = find(Keep(:,i));
    indl = nlinks_e(i);    
    if indl == 0, continue; end    
    k = k + indl;    
    events.event(j:k) = i;    
    events.trajIndex(j:k) = ind;
    j = k + 1;    
end

events = join(events, indexChange);
events.trajIndex = [];
events.Properties.VariableNames{2} = 'trajIndex';


% Filter the forcing data
for i = 1:nfields
    s = size(Forc.(fields{i}));
    sl = length(s);
    y = s == nTraj;
    if any(y)
        if sl == 1
            Forc.(fields{i}) = Forc.(fields{i})(iTraj);
        elseif sl == 2
            Forc.(fields{i}) = Forc.(fields{i})(:,iTraj);
        elseif sl == 3
            Forc.(fields{i}) = Forc.(fields{i})(:,:,iTraj);
        end
    end
end




