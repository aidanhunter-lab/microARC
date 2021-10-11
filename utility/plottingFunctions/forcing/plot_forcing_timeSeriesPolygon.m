function plot_forcing_timeSeriesPolygon(var, fixedParams, forcing, dat, varargin)

extractVarargin(varargin)

if ~exist('depth', 'var')
    depth = 'mean';
end

if ~exist('polygonColour', 'var')
    polygonColour = [0 0 0.5];
end
if ~exist('lineColour', 'var')
    lineColour = [0 0 0];
end

if ~exist('highlightColour', 'var')
    highlightColour = [0 1 0]; % green as default
end
if ~exist('plotNew', 'var')
    plotNew = true;
end
if ~exist('fixedYaxis', 'var')
    fixedYaxis = false;
end

if ~exist('axesTextSize', 'var')
    axesTextSize = 12;
end
if ~exist('titleTextSize', 'var')
    titleTextSize = 10;
end
if ~exist('legendTextSize', 'var')
    legendTextSize = 10;
end
if ~exist('legendPointSize', 'var')
    legendPointSize = 36;
end

if exist('waterMass', 'var')
    waterMass = eval('waterMass');
    if ~any(strcmp(waterMass, forcing.waterMass))
        error('Optional argument "waterMass" does not match values in forcing data.')
    end
    traj = find(strcmp(waterMass, forcing.waterMass));
else
%     waterMass = [];
    if ~exist('traj', 'var')
        traj = 1:forcing.nTraj;
    end
end
singleTrajectory = length(traj) == 1;

if singleTrajectory
    warning('Only a single trajectory has been selected. Multiple trajectories are required to define polygon.')
end

%%

x = yearday(forcing.t(:,1)); % yearday on x-axis
% if exist('event', 'var')
%     event = eval('event');
%     etime = unique(dat.scalar.Yearday(dat.scalar.Event == event,:)); % yearday of sampling event
% end


switch var
    case 'temperature'
        y = forcing.T(:,:,traj);
        ylab = ['temperature (' char(176) 'C)'];
    case 'diffusivity'
        y = forcing.K_center(:,:,traj);
        ylab = 'diffusivity (m^2 day^{-1})';
    case 'PAR'
        y = forcing.PARsurf(:,:,traj);
        y = 1e-6 .* y;
        ylab = 'PAR (Ein day^{-1} m^{-2})';
end

y = depthAggregate(y, depth, fixedParams);

switch singleTrajectory
    case true % line plot
        plot(x, y, 'Color', [0 0 0])
    case false % polygon
        lo = min(y, [], ndims(y)); % min and max over trajectories
        hi = max(y, [], ndims(y));
        ym = median(y, ndims(y)); % median over trajectories
        lo = lo(:); hi = hi(:); ym = ym(:);
        tt = unique(table([x; flip(x)], [lo; flip(hi)]), 'stable');
        pgon = polyshape(tt{:,1}, tt{:,2});
        plot(pgon, 'FaceColor', polygonColour);
        switch var
            case 'temperature', YLim = [min(lo) max(hi)];
            otherwise, YLim = [0 max(hi)];
        end
        set(gca, {'XLim', 'YLim'}, {[min(x) max(x)] , YLim})
        hold on
        plot(x,ym,'Color',lineColour)
        hold off
end

if ~exist('xlab', 'var')
    xlab = 'year-day';
end
if strcmp(depth, 'surface')
    ylab = {ylab, 'surface layer'};
elseif strcmp(depth, 'mean')
    ylab = {ylab, 'depth averaged'};
elseif strcmp(depth, 'max')
    ylab = {ylab, 'max over depths'};
end
xlabel(xlab)
ylabel(ylab)



end


function v = depthAggregate(x, depth, fixedParams)
xsize = size(x);
if strcmp(depth, 'surface')
    v = reshape(x(1,:), xsize(2:end));
end
if strcmp(depth, 'mean')
    v = squeeze(sum(fixedParams.zwidth .* x) ./ fixedParams.Htot);
end
if strcmp(depth, 'max')
    v = squeeze(max(x));
end
end
