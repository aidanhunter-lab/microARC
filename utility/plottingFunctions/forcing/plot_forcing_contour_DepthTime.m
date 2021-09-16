function plot_forcing_contour_DepthTime(var, forcing, auxVars, ...
    fixedParams, varargin)

extractVarargin(varargin)

if ~exist('traj', 'var')
    traj = 1:forcing.nTraj; % By default include all trajectories in plot after averaging
end

if ~exist('logScale', 'var')
    logScale = false; % log scale for plotted variable?
end

if ~exist('scaleFactor', 'var')
    scaleFactor = 1;
end

if ~exist('smooth', 'var')
    smooth = 'linear'; % interpolation method for contour plot
end

if ~exist('ensurePositive', 'var')
    switch var
        case 'T', ensurePositive = false;
        otherwise, ensurePositive = true; % interpolation can produce negative values which are reset
    end
end

if ~exist('maxDepth', 'var')
    maxDepth = abs(fixedParams.zw(end)); % maximum depth to plot - if unspecified then defaults to maximum modelled depth
end

if ~exist('ColourBar', 'var')
    ColourBar = true;
end

if ~exist('ColourBarLabel', 'var')
    ColourBarLabel = '';
end

if ~exist('yLabel', 'var')
    yLabel = 'depth (m)';
end
if ~exist('xLabel', 'var')
    xLabel = 'time (yearday)';
end

if ~exist('nYTicks', 'var')
    nYTicks = 7; % number of ticks on y-axis
end

if ~exist('TitleWeight', 'var')
    TitleWeight = 'normal';
end

% Time-depth grid for interpolation
[depth, time] = ndgrid(abs(fixedParams.z), 1:fixedParams.nt);
[depthGrid, timeGrid] = ndgrid(1:1:abs(fixedParams.zw(end)), 1:fixedParams.nt);
depthGrid = depthGrid(1:maxDepth,:);
timeGrid = timeGrid(1:maxDepth,:);

switch var
    case 'T'
        x = forcing.T(:,:,traj);
    case 'K'
        x = forcing.K_center(:,:,traj);
    case {'PAR','I'}
        x = squeeze(auxVars.I(:,:,:,traj));
end

x = scaleFactor .* x;

ntraj = length(traj);
if ntraj > 1
    x = mean(x, ndims(x), 'omitnan'); % some trajectories may include NaNs from encountering shallow seafloors -- ignore these values when averaging over multiple trajectories
end
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
if ensurePositive, Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0)); end

switch logScale, case true
    Fsmooth = log10(Fsmooth);
end

contourf(Fsmooth)

if ColourBar
    cb = colorbar;
    cb.Label.String = ColourBarLabel;
    switch logScale, case true
        for ii = 1:length(cb.TickLabels), cb.TickLabels{ii} = string(round(10 ^ str2double(cb.TickLabels{ii}),2,'significant')); end
    end
end

if exist('Title', 'var')
    title(Title, 'FontWeight', TitleWeight)
end

xlabel(xLabel)
ylabel(yLabel)

xticks(100:100:fixedParams.nt)
xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))

yt = round(linspace(0,max(abs(depthGrid(:))),nYTicks), 2, 'significant');
yticks(yt)
% yt(2:end) = -yt(2:end); % display depths as negative
yticklabels(flip(yt))




% switch var
%     
%     case 'forcing'
%         
%         plt = figure;
%         plt.Units = 'inches';
%         plt.Position = [0 0 8 9];
%         
%         % Temperature
%         subplot(3,1,1)
%         x = forcing.T(:,:,traj);
%         if ntraj > 1
%             x = mean(x, ndims(x));
%         end
%         F = griddedInterpolant(depth, time, x, smooth);
%         
%         Fsmooth = flip(F(depthGrid, timeGrid));
% %         Fsmooth = flip(F(depthGrid, timeGrid));
%         contourf(Fsmooth)
%         cb = colorbar;
%         cb.Label.String = '\circC';
%         title('Temperature')
%         ylabel('depth (m)')
%         xticks(100:100:fixedParams.nt)
%         xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
%         
%         yt = round(linspace(0,max(abs(depthGrid(:))),7), 2, 'significant');
%         yticks(yt)
%         yticklabels(flip(-yt))
%         
%         
%         % Diffusivity
%         subplot(3,1,2);
%         x = forcing.K_center(:,:,traj);
%         if ntraj > 1
%             x = mean(x, ndims(x));
%         end
%         F = griddedInterpolant(depth, time, x, smooth);
%         Fsmooth = flip(F(depthGrid, timeGrid));
%         Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0));
% %         Fsmooth(Fsmooth<=0) = nan;
%         contourf(log10(Fsmooth))
%         cb = colorbar;
%         for ii = 1:length(cb.TickLabels), cb.TickLabels{ii} = string(round(10 ^ str2double(cb.TickLabels{ii}),2,'significant')); end
%         cb.Label.String = 'm^2 day^{-1}';
%         title('Diffusivity')
%         ylabel('depth (m)')
%         xticks(100:100:fixedParams.nt)
%         xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))
%         
%         yt = round(linspace(0,max(abs(depthGrid(:))),7), 2, 'significant');
%         yticks(yt)
%         yticklabels(flip(-yt))
%         %                 set(gca,'colorscale','log')
%         
%         % PAR
%         subplot(3,1,3)
%         x = squeeze(auxVars.I(:,:,:,traj));
%         if ntraj > 1
%             x = mean(x, ndims(x), 'omitnan');
%         end
%         F = griddedInterpolant(depth, time, x, smooth);
%         Fsmooth = flip(F(depthGrid, timeGrid));
%         Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0));
% %         Fsmooth(Fsmooth<=0) = nan;
%         contourf(Fsmooth)
%         cb = colorbar;
%         cb.Label.String = '\muEin day^{-1} m^{-2}';
%         title('PAR')
%         xlabel('year-day')
%         ylabel('depth (m)')
%         xticks(100:100:fixedParams.nt)
%         xticklabels(yearday(forcing.t(100:100:fixedParams.nt,1)))        
%         yt = round(linspace(0,max(abs(depthGrid(:))),7), 2, 'significant');
%         yticks(yt)
%         yticklabels(flip(-yt))
% %         yticks(linspace(0,abs(fixedParams.zw(end)),7))
% %         yticklabels(linspace(fixedParams.zw(end),0,7))
%         colormap plasma
%         if exist('sampleEvent', 'var')
%             sgtitle(['Sample event ' num2str(sampleEvent) ': ' waterOrigin ' origin'])
%         else
%             sgtitle([waterOrigin ' origin'])
%         end
%         
        %----------------------------------------------------------
        
        