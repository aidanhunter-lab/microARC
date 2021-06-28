function [Forc, events] = chooseTrajectories(Forc, dat, maxDist, numTraj)

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
    sLat = dat.Latitude(ie);
    sLon = dat.Longitude(ie);
    sTime = dat.Yearday(ie);
    samplePos_rads = deg2rad([sLon sLat]);

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
    
    Keep(:,i) = keep;

end

% For each sampling event choose numTraj particle trajectories from those 
% within maxDist of the event.
% Use a clustering method to choose the most dissimilar trajectories,
% thereby maximising variability in the forcing data.

for i = 1:nEvent
    i_traj = Keep(:,i);
    f_traj = find(i_traj);
    n_traj = sum(i_traj);
    r_traj = Rank(f_traj,i);
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
        clust = cluster(lat_link,'maxclust',numTraj);
        
        % From each cluster select the trajectory closest to sample event i.
        % If there are too few trajectories within radius maxDist of sample
        % location then include duplicates to return numTraj trajectories.
        uc = unique(clust);
        nc = length(uc);
        useTrajectories = nan(numTraj,1);
        for j = 1:nc
            jc = clust == j;
            trajClust = f_traj(jc); % trajectories in cluster j
            trajRank = r_traj(jc);
            useTrajectories(j) = trajClust(trajRank == min(trajRank)); % from cluster j choose trajectory closest to sampling event i
        end
        needMore = any(isnan(useTrajectories));
        if needMore
            warning(['Duplicate trajectories used for sampling event ' ...
                num2str(i) '. This is nothing to worry about unless warning' ...
                ' is repeated for lots of event numbers, in which case try' ...
                ' increasing maxDist or decreasing numTraj.'])
        end
        ind = 1;
        while needMore % include duplicates if needed
            useTrajectories(find(isnan(useTrajectories), 1)) = ... 
                f_traj(r_traj == ind);
            ind = ind + 1;
            if ind > max(r_traj), ind = 1; end
            needMore = any(isnan(useTrajectories));
        end
        Keep(:,i) = zeros(size(Keep,1), 1);
        for j = 1:numTraj
            Keep(useTrajectories(j),i) = Keep(useTrajectories(j),i) + 1;
        end
    else
        Keep(:,i) = zeros(size(Keep,1), 1);
        Keep(f_traj,i) = 1;
    end
end


iTraj = any(Keep,2); % index selected trajectories
nTraj_new = sum(iTraj); % number of unique trajectories selected
Forc.nTraj = nTraj_new;
fields = fieldnames(Forc);
nfields = length(fields);
% Forc = orderfields(Forc, [1 nfields 2:nfields-1]);

trajIndex = find(iTraj);
trajIndex_new = 1:nTraj_new;
trajIndex = trajIndex(:); trajIndex_new = trajIndex_new(:);
indexChange = table(trajIndex, trajIndex_new); % filtering trajectories changes the indexing

% store event numbers and associated trajectories in a table
nlinks_e = sum(Keep);
if sum(nlinks_e == 0) > 1
    warning(['Sampling events ' num2str(find(nlinks_e == 0)) ...
        ' were discarded as no trajectories were within a distance of' ...
        ' maxDist from sample location. Increasing maxDist might help,' ...
        ' although some sampling event locations are far from any' ...
        ' trajectories...'])
elseif sum(nlinks_e == 0) == 1
    warning(['Sampling event ' num2str(find(nlinks_e == 0)) ...
        ' was discarded as no trajectories were within a distance of' ...
        ' maxDist from sample location. Increasing maxDist might help,' ...
        ' although some sampling event locations are far from any' ...
        ' trajectories...'])
end


nlinks = sum(nlinks_e);
events = table(nan(nlinks,1),nan(nlinks,1));
events.Properties.VariableNames = {'event', 'trajIndex'};
j=1; k=0;
for i = 1:nEvent    
    ind = find(Keep(:,i));
    indrep = Keep(ind,i); % repeat for any duplicated trajectories
    ind = arrayfun(@(x) repmat(ind(x), indrep(x), 1), 1:numel(ind), 'uni', 0);
    ind = vertcat(ind{:});
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

ue = unique(events.event);
ne = length(ue);
for i = 1:ne % catch any events with a single trajectory
    e = ue(i);
    ie = events.event == e;
    if sum(ie) == 1
        p = events(ie,:);
        events(ie,:) = [];
        events = [events; repmat(p, [numTraj, 1])];
    end
end
[~,o] = sort(events.event);
events = events(o,:);


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
