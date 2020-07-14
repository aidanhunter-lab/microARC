function [Forc, events] = chooseTrajectories(Forc, dat, maxDist, maxTraj)

% Select trajectories for each sampling event based on distance between
% particles and sample.
% Filter out any unused forcing data.
% Return a table linking event numbers to particle indexes.

% nt = size(Forc.t,1);
nTraj = length(Forc.iTraj);

% particle times as year and yearday
time = yearday(Forc.t);
% [year,~] = datevec(Forc.t);

% For each event, filter out all particles further than maxDist from event
% location at event time.
nEvent = length(unique(dat.Event));
Dist = nan(nTraj, nEvent);
Rank = nan(size(Dist));
Keep = nan(size(Dist));

for i = 1:nEvent
    ie = find(dat.Event == i,1);
%     sYear = dat.Year(ie);
    sLat = dat.Latitude(ie);
    sLon = dat.Longitude(ie);
    sTime = dat.Yearday(ie);
    samplePos_rads = deg2rad([sLon sLat]);

%     keep2d = sTime == time & sYear == year;    
    keep2d = sTime == time;    
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
    
%     % Limit number of trajectories to maxTraj per event
%     nbest = Rank(:,i) <= maxTraj;
%     keep = keep(:) & nbest(:);
%     Dist(~keep,i) = nan;
%     Rank(~keep,i) = nan;
    
    Keep(:,i) = keep;
    
end

% For each sampling event choose, at most, maxTraj particle trajectories
% from those within maxDist of the event.
% Use a clustering method to choose the most dissimilar trajectories,
% thereby maximising variability in the forcing data.

for i = 1:nEvent
    i_traj = Keep(:,i);
    f_traj = find(i_traj);
    n_traj = sum(i_traj);
    if n_traj == 0, continue; end    
    if n_traj ~= 1
        % We're only interested in measuring trajectory dissimilarities
        % leading up to the sampling event... where the trajectories travel
        % after that is irrelevant.
        times = 1:unique(dat.Yearday(dat.Event == i));
        dist_dtw = zeros(n_traj); % distances between particles
        for j = 2:n_traj
            for k = 1:j-1
                s1 = [Forc.T(1,times,f_traj(j))' Forc.K(1,times,f_traj(j))' ...
                    Forc.PARsurf(1,times,f_traj(j))'];
                s2 = [Forc.T(1,times,f_traj(k))' Forc.K(1,times,f_traj(k))' ...
                    Forc.PARsurf(1,times,f_traj(k))'];
                dist_dtw(j,k) = dtw(s1,s2);
            end
        end
        lat_link = linkage(dist_dtw);
        clust = cluster(lat_link,'maxclust',maxTraj);
        % Select the 1st trajectory from each cluster
        uc = unique(clust);
        nc = length(uc);
        useTrajectories = nan(nc,1);
        for j = 1:nc
            useTrajectories(j) = f_traj(find(clust == j, 1));
        end
        Keep(:,i) = false;
        Keep(useTrajectories,i) = true;
    else
        Keep(:,i) = false;
        Keep(f_traj,i) = true;
    end
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




