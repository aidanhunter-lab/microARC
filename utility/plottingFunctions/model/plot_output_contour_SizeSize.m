function plot_output_contour_SizeSize(var, Out, AuxVars, ...
    fixedParams, forcing, varargin)

% Contour plots of predation losses integrated over depth and time

if forcing.nTraj ~= size(Out.N, ndims(Out.N))
    error('Number of trajectories contained in forcing data struct must match the number of trajectories contained in model output struct.')
end

extractVarargin(varargin)

if ~exist('xlab', 'var')
    xlab = 'ESD (\mum)';
end
if ~exist('ylab', 'var')
    ylab = 'ESD (\mum)';
end
if ~exist('smooth', 'var')
    smooth = 'linear';
end
if ~exist('ensurePositive', 'var')
    ensurePositive = true; % interpolation can produce negative values which are reset
end

if ~exist('Times', 'var')
    Times = 1:fixedParams.nt;
end

if ~exist('totOrAve', 'var')
    totOrAve = 'total';
end

if ~exist('logScale', 'var')
    logScale = false;
else
    logScale = eval('logScale');
end
if ~exist('powerScale', 'var')
    powerScale = [];
else
    powerScale = eval('powerScale');
end

% Index trajectories to plot
if exist('waterOrigin', 'var')
    waterOrigin = eval('waterOrigin');
    if ~contains(waterOrigin, forcing.waterMass)
        error('Optional argument "waterOrigin" must match one of the "waterMass" values in forcing data.')
    end
    traj = find(strcmp(forcing.waterMass, waterOrigin));
elseif exist('traj', 'var')
    traj = eval('traj');
    if ~islogical(traj)
        if any(~ismember(traj, 1:forcing.nTraj))
            error('Values in optional indexing argument "traj" must not exceed number of trajectories stored in forcing data.')
        end
    else
        if length(traj) ~= forcing.nTraj
            error('Length of logical indexing argument "traj" must equal number of trajectories in forcing data.')
        end
        traj = find(traj);
    end
else
    traj = 1:forcing.nTraj;
end


% Extract outputs
vars = {'Z_P_C', 'Z_P_N', 'Z_Z_C', 'Z_Z_N'};
nvars = length(vars);

losses_depthInt = squeeze(AuxVars.predation_losses_depthInt);
losses_depthInt = losses_depthInt(:,:,:,Times,traj);

dt = diff(forcing.t(1:2,1));
nt = length(Times);

losses_ave = squeeze(mean(losses_depthInt, find(size(losses_depthInt) == nt), 'omitnan'));
losses_tot = squeeze(sum(dt .* losses_depthInt, find(size(losses_depthInt) == nt), 'omitnan'));

switch totOrAve
    case 'total', losses = losses_tot;
    case 'average', losses = losses_ave;
end

% losses_tot = squeeze(AuxVars.predation_losses_tot);
% losses_tot = losses_tot(:,:,:,traj);

Z_P_C = squeeze(losses(:,1:fixedParams.nPP_size,fixedParams.ZP_C_index,:));
Z_P_N = squeeze(losses(:,1:fixedParams.nPP_size,fixedParams.ZP_N_index,:));
Z_Z_C = squeeze(losses(:,fixedParams.nPP_size+1:end,fixedParams.ZP_C_index,:));
Z_Z_N = squeeze(losses(:,fixedParams.nPP_size+1:end,fixedParams.ZP_N_index,:));

x = nan([size(Z_P_C) nvars]);
for i = 1:nvars
    x(:,:,:,i) = eval(vars{i});
end

% If multiple trajectories selected then take averages
ntraj = length(traj);
if ntraj > 1
    x = squeeze(mean(x, ndims(x)-1, 'omitnan'));
end

x = x(:,:,strcmp(var, vars));

if ~isempty(powerScale)
    x = x .^ powerScale;
end

if logScale
   x = log10(x+1);
%    ensurePositive = false;
end

pred = var(1);
prey = var(3);

% Time-size grid for interpolation
[SizePred, SizePrey] = ndgrid(abs(fixedParams.([pred 'Pdia'])), abs(fixedParams.([prey 'Pdia'])));
[predGrid, preyGrid] = ndgrid(1:1:fixedParams.([pred 'Pdia_intervals'])(end), 1:1:fixedParams.([prey 'Pdia_intervals'])(end));

if exist('xLim', 'var')
    predGrid = predGrid(:,xLim(1):xLim(2));
    preyGrid = preyGrid(:,xLim(1):xLim(2));
end
if exist('yLim', 'var')
    predGrid = predGrid(yLim(1):yLim(2),:);
    preyGrid = preyGrid(yLim(1):yLim(2),:);
end

F = griddedInterpolant(SizePred, SizePrey, x, smooth);
Fsmooth = F(predGrid, preyGrid);
if ensurePositive, Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0)); end

if exist('levels', 'var')
    contourf(Fsmooth, levels)
else
   contourf(Fsmooth) 
end

set(gca, {'XScale', 'YScale'}, {'log', 'log'})

cb = colorbar;

if ~isempty(powerScale)
    TickLabels = get(cb, 'TickLabels');
    TickLabels = cellfun(@(z) str2num(z), TickLabels) .^ (1 / powerScale);
    TickLabels = round(TickLabels, 2 , 'significant');
    TickLabels = arrayfun(@(z) num2str(z), TickLabels, 'UniformOutput', false);
    set(cb, 'TickLabels', TickLabels)
end

% if logScale
%     if ~exist('colBarMax', 'var')
%         colBarMax = 0.5;
%     end
%     Ticks = get(cb, 'Ticks');
%     Ticks = unique([Ticks, Ticks(end) + colBarMax*(max(x(:))-Ticks(end))]);
% %     Ticks = [Ticks, max(x(:))];
%     Ticks = round(10 .^ Ticks, 2, 'significant');
%     TickLabels = arrayfun(@(z) num2str(z), Ticks(:), 'UniformOutput', false);
%     Ticks = log10(Ticks);
%     set(cb, {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
% end

if logScale
    Ticks = get(cb, 'Ticks');
    tl = round(-1 + 10 .^ Ticks([1,end]), 1, 'significant');
    ltl = log10(tl+1);
    Ticks = linspace(ltl(1), ltl(2), length(Ticks));
    TickLabels = round(-1 + 10 .^ Ticks, 2, 'significant');
    TickLabels = arrayfun(@(z) num2str(z), TickLabels(:), 'UniformOutput', false);
    set(cb, {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
end

if exist('ColourBarLabel','var')
   cb.Label.String = eval('ColourBarLabel') ;
end

xlabel(xlab)
ylabel(ylab)
if ~exist('Title', 'var')
    Title = [];
end
title(Title)

if exist('xLim', 'var')
    XTick = get(gca, 'XTick');
    set(gca, 'XTickLabel', arrayfun(@(z) num2str(z), XTick + xLim(1), 'UniformOutput', false))
end



        