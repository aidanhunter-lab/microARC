function [Forc, Data] = particleOrigin(forc, data, varargin)

% Do data values depend on whether samples are from Arctic or Atlantic waters?
% Associate each sample event with either Arctic or Atlantic water.
% Use dynamic time warping metrics to create matrix of distances separating
% forcing data particle trajectories, then use hierchical clustering to
% pick out the 2 main clusters which will correspond to Arctic/Atlantic

Forc = forc;
Data = data;

% extra arguments can be name-value pairs controlling whether to display
% plots
v = reshape(varargin, [2 0.5*length(varargin)]);

dist_dtw = zeros(Forc.nTraj);
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
end
lat_link = linkage(dist_dtw);

if ~isempty(v) && any(strcmp(v(1,:),'dendrogram')) && v{2,strcmp(v(1,:),'dendrogram')}
    figure
    dendrogram(lat_link,'ColorThreshold','default')
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

% Plot result of clustering
if ~isempty(v) && any(strcmp(v(1,:),'trajectoryPlot')) && v{2,strcmp(v(1,:),'trajectoryPlot')}
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
waterMass = zeros(Data.nEvents,2); % for each event, count particles originating from Atlantic or Arctic
for i = 1:Data.nEvents
    p = Data.EventTraj(i,:); % particles used for event i
    x = Forc.waterMass(p);
    waterMass(i,1) = sum(strcmp(x,'Atlantic'));
    waterMass(i,2) = sum(strcmp(x,'Arctic'));
end
% Given the way in which trajectories have been filtered simply by distance
% from sampling event, it is likely that most (or all) trajectories used
% for each event entirely originate from either the Arctic or Atlantic, so
% let's just assume there is zero mixing so that each sampling event is
% either within Atlantic OR Arctic water... This could/should be modified
% later...
ind = waterMass ./ sum(waterMass,2) > 0.5;
Data.waterMass = cell(Data.nEvents,1);
Data.waterMass(ind(:,1)) = {'Atlantic'};
Data.waterMass(ind(:,2)) = {'Arctic'};

