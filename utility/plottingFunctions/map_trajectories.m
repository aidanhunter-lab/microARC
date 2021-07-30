% create map of  trajectory plots, optionally highlight some trajectories

function plt = map_trajectories(forc, fixedParams, varargin)

% extract trajectories
latTraj = forc.y;
lonTraj = forc.x; 
timeTraj = forc.t; 
idTraj = forc.iTraj; 

% plot
plt = figure
for i = 1:forc.nTraj
    if i == 2 
        hold on 
    end
    
    latTraj = forc.y(:,i);
    lonTraj = forc.x(:,i); 
    timeTraj = yearday(forc.t(:,i)); 
    idTraj = repelem(forc.iTraj(i), fixedParams.nt);
    waterMass = forc.waterMass(i); 
    
    if strcmp(waterMass,'Atlantic')
        col = [1 0 0 0.1];
    else
        col = [0 0 1 0.1];
    end
    
    g = geoplot(latTraj, lonTraj, '-', 'Color', col, 'lineWidth', 0.5);
    g.DataTipTemplate.DataTipRows(3) = dataTipTextRow('ID', idTraj);
    g.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Yearday', timeTraj);
end

% optional params:
% highlight specific trajectories?
if ~isempty(varargin)
    if any(strcmp(varargin, 'highlightTrajs'))
        highlightTrajs = varargin{find(strcmp(varargin, 'highlightTrajs'))+1};
     
        % print trajectories to be highlighted
        highlightTrajs
        
        % find indices of these trajectories 
        trajs = find(ismember(forc.iTraj, highlightTrajs));
        
        % add them to plot
        hold on
        for i = trajs
            latTraj = forc.y(:,i);
            lonTraj = forc.x(:,i); 
            timeTraj = yearday(forc.t(:,i));
            idTraj = repelem(forc.iTraj(i), fixedParams.nt);
            waterMass = forc.waterMass(i); 
            if strcmp(waterMass,'Atlantic')
                col = 'r';
            else
                col = 'b';
            end
            
            g = geoplot(latTraj, lonTraj, '-', 'Color', col, 'lineWidth', 0.5) % , ...
            g.DataTipTemplate.DataTipRows(3) = dataTipTextRow('ID', idTraj);
            g.DataTipTemplate.DataTipRows(4) = dataTipTextRow('Yearday', timeTraj);
        end
        hold off
    end
end



end

