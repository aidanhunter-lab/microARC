function [Forc, Data] = particleOrigin(Forc, varargin)
% Do data values depend on whether samples are from Arctic or Atlantic waters?
% Associate each sample event with either Arctic or Atlantic water.
% Use dynamic time warping metrics to create matrix of distances separating
% forcing data particle trajectories, then use hierchical clustering to
% pick out the 2 main clusters which will correspond to Arctic/Atlantic

extractVarargin(varargin)

if ~exist('trajectoryPlot', 'var')
    trajectoryPlot = false;
end

if ~exist('dendrogramPlot', 'var')
    dendrogramPlot = false;
end

if ~exist('progressBar', 'var')
    progressBar = false;
end

Forc.nTraj = length(Forc.iTraj);

calculateGroups = ~isfield(Forc, 'waterMass') || dendrogramPlot;

if calculateGroups
    
    dist_dtw = zeros(Forc.nTraj);
    switch progressBar, case true, progress = waitbar(0, 'Group forcing data by particle origin'); end
    for i = 2:Forc.nTraj
        for j = 1:i-1
            % Use particle latitudes to separate Atlantic and Arctic waters
            s1 = Forc.y(:,i);
            s2 = Forc.y(:,j);
            dist_dtw(i,j) = dtw(s1,s2);
            if i == Forc.nTraj && j == i-1
                dist_dtw = dist_dtw + dist_dtw';
            end
        end
        switch progressBar, case true, waitbar((i-1)/(Forc.nTraj-1) , progress); end
    end
    switch progressBar, case true, delete(progress); end
    lat_link = linkage(dist_dtw);
    
    if exist('dendrogramPlot', 'var') && dendrogramPlot
        figure
        dendrogram(lat_link,'ColorThreshold','default');
        xlabel('particle groups')
        title({'Particle trajectories clustered by latitude', ...
            'Two main clusters separate water of Arctic & Atlantic origin'})
    end
    
    clust = cluster(lat_link,'maxclust',2);
    waterMass = cell(size(clust));
    if mean(Forc.y(1,clust==1)) < mean(Forc.y(1,clust==2))
        waterMass(clust==1) = {'Atlantic'}; waterMass(clust==2) = {'Arctic'};
    else
        waterMass(clust==1) = {'Arctic'}; waterMass(clust==2) = {'Atlantic'};
    end
    
    Forc.waterMass = waterMass'; % Does particle originate from Arctic or Atlantic ocean?
    
end

% Plot result of clustering
if exist('trajectoryPlot', 'var') && trajectoryPlot
    figure
    for i = 1:Forc.nTraj
        if i==2
            hold on
        end
        if strcmp(Forc.waterMass{i}, 'Atlantic')
            plot(Forc.y(:,i),'r')
        else
            plot(Forc.y(:,i),'b')
        end
        if i==Forc.nTraj
            xlabel('year-day')
            ylabel(['latitude (' char(176) 'N)'])
            title('Particle trajectories clustered by latitude')
            x = gca;
            xl = x.XLim;
            yl = x.YLim;
            text(xl(1) + 3/5 * diff(xl), yl(1) + 0.95 * diff(yl), 'Arctic water')
            plot(xl(1) + [0.5 11/20] * diff(xl), repmat(yl(1) + 0.95 * diff(yl), [1 2]), 'b')
            text(xl(1) + 3/5 * diff(xl), yl(1) + 0.90 * diff(yl), 'Atlantic water')
            plot(xl(1) + [0.5 11/20] * diff(xl), repmat(yl(1) + 0.90 * diff(yl), [1 2]), 'r')
            hold off
        end
    end
end

% Associate each sample with Arctic or Atlantic waters, or both...
if exist('Data', 'var')
    Data = eval('Data');    
    waterMass = zeros(Data.scalar.nEvents,2); % for each event, count particles originating from Atlantic or Arctic
    for i = 1:Data.scalar.nEvents
        p = Data.scalar.EventTraj(i,:) > 0; % particles used for event i
        x = Forc.waterMass(p);
        waterMass(i,1) = sum(strcmp(x,'Atlantic'));
        waterMass(i,2) = sum(strcmp(x,'Arctic'));
    end
    Data.scalar.waterMass = cell(Data.scalar.nEvents,1); % label each sampling event as from Arctic, Atlantic, or a mix of water depending on origin of particle trajectories
    Arctic = waterMass(:,1) == 0;
    Atlantic = waterMass(:,2) == 0;
    Mix = ~(Arctic | Atlantic);
    Data.scalar.waterMass(Arctic) = {'Arctic'};
    Data.scalar.waterMass(Atlantic) = {'Atlantic'};
    Data.scalar.waterMass(Mix) = {'Arctic/Atlantic'};
    Data.scalar.AtlanticOrigin = waterMass(:,1) ./ sum(waterMass,2); % proportion of trajectories of Atlantic or Arctic origin for each sampling event
    Data.scalar.ArcticOrigin = waterMass(:,2) ./ sum(waterMass,2);
    
    % same for size data
    events = unique(Data.sizeFull.Event);
    nEvents = length(events);
    waterMass = zeros(nEvents,2);
    for i = 1:nEvents
        eventi = events(i);
        p = Data.scalar.EventTraj(eventi,:) > 0;
        x = Forc.waterMass(p);
        waterMass(i,1) = sum(strcmp(x,'Atlantic'));
        waterMass(i,2) = sum(strcmp(x,'Arctic'));
    end
    
    Data.sizeFull.waterMass = cell(nEvents,1);
    Arctic = waterMass(:,1) == 0;
    Atlantic = waterMass(:,2) == 0;
    Mix = ~(Arctic | Atlantic);
    Data.sizeFull.waterMass(Arctic) = {'Arctic'};
    Data.sizeFull.waterMass(Atlantic) = {'Atlantic'};
    Data.sizeFull.waterMass(Mix) = {'Arctic/Atlantic'};
    Data.sizeFull.AtlanticOrigin = waterMass(:,1) ./ sum(waterMass,2); % proportion of trajectories of Atlantic or Arctic origin for each sampling event
    Data.sizeFull.ArcticOrigin = waterMass(:,2) ./ sum(waterMass,2);
    
    Data.sizeFull.dataBinned.waterMass = Data.sizeFull.waterMass;
    Data.sizeFull.dataBinned.AtlanticOrigin = Data.sizeFull.AtlanticOrigin;
    Data.sizeFull.dataBinned.ArcticOrigin = Data.sizeFull.ArcticOrigin;
    
end
