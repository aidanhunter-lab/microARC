function plot_output_contour_SizeTime(var, Out, AuxVars, ...
    fixedParams, forcing, varargin)

if forcing.nTraj ~= size(Out.N, ndims(Out.N))
    error('Number of trajectories contained in forcing data struct must match the number of trajectories contained in model output struct.')
end

extractVarargin(varargin)

if ~exist('xlab', 'var')
    xlab = 'year-day';
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

if ~exist('logScale', 'var')
    logScale = false;
else
    logScale = eval('logScale');
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
vars = {'P_N', 'P_Chl', 'P_C', 'Z_N', 'Z_C', 'P_cellDensity', 'Z_cellDensity', 'P_biovolume', 'Z_biovolume', ...
    'P_ON', 'P_OC', 'Z_ON', 'Z_OC', 'P_DOC', 'P_POC', 'P_DON', 'P_PON', ...
    'Z_DOC_mort', 'Z_POC_mort', 'Z_DOC_mess', 'Z_POC_mess', 'Z_DON_mort', 'Z_PON_mort', 'Z_DON_mess', 'Z_PON_mess'};
nvars = length(vars);
P_C = squeeze(Out.P(:,:,fixedParams.PP_C_index,:,traj));
P_N = squeeze(Out.P(:,:,fixedParams.PP_N_index,:,traj));
P_Chl = squeeze(Out.P(:,:,fixedParams.PP_Chl_index,:,traj));
Z_C = squeeze(Out.Z(:,:,fixedParams.ZP_C_index,:,traj));
Z_N = squeeze(Out.Z(:,:,fixedParams.ZP_N_index,:,traj));
P_cellDensity = AuxVars.cellDensity(1:fixedParams.nPP_size,:,:,traj);
Z_cellDensity = AuxVars.cellDensity(fixedParams.nPP_size+1:end,:,:,traj);
P_biovolume = AuxVars.biovolume(1:fixedParams.nPP_size,:,:,traj);
Z_biovolume = AuxVars.biovolume(fixedParams.nPP_size+1:end,:,:,traj);

P_DON_mort = squeeze(AuxVars.OM_mort_DOM(1:fixedParams.nPP_size,:,fixedParams.PP_N_index,:,traj));
P_DOC_mort = squeeze(AuxVars.OM_mort_DOM(1:fixedParams.nPP_size,:,fixedParams.PP_C_index,:,traj));
Z_DON_mort = squeeze(AuxVars.OM_mort_DOM(fixedParams.nPP_size+1:end,:,fixedParams.PP_N_index,:,traj));
Z_DOC_mort = squeeze(AuxVars.OM_mort_DOM(fixedParams.nPP_size+1:end,:,fixedParams.PP_C_index,:,traj));
Z_DON_mess = squeeze(AuxVars.OM_mess_DOM(:,:,:,fixedParams.PP_N_index,:,traj));
Z_DOC_mess = squeeze(AuxVars.OM_mess_DOM(:,:,:,fixedParams.PP_C_index,:,traj));
P_PON_mort = squeeze(AuxVars.OM_mort_POM(1:fixedParams.nPP_size,:,fixedParams.PP_N_index,:,traj));
P_POC_mort = squeeze(AuxVars.OM_mort_POM(1:fixedParams.nPP_size,:,fixedParams.PP_C_index,:,traj));
Z_PON_mort = squeeze(AuxVars.OM_mort_POM(fixedParams.nPP_size+1:end,:,fixedParams.PP_N_index,:,traj));
Z_POC_mort = squeeze(AuxVars.OM_mort_POM(fixedParams.nPP_size+1:end,:,fixedParams.PP_C_index,:,traj));
Z_PON_mess = squeeze(AuxVars.OM_mess_POM(:,:,:,fixedParams.PP_N_index,:,traj));
Z_POC_mess = squeeze(AuxVars.OM_mess_POM(:,:,:,fixedParams.PP_C_index,:,traj));
P_DON = P_DON_mort;
P_DOC = P_DOC_mort;
P_PON = P_PON_mort;
P_POC = P_POC_mort;
Z_DON = Z_DON_mort+Z_DON_mess;
Z_DOC = Z_DOC_mort+Z_DOC_mess;
Z_PON = Z_PON_mort+Z_PON_mess;
Z_POC = Z_POC_mort+Z_POC_mess;
P_ON = P_DON + P_PON;
P_OC = P_DOC + P_POC;
Z_ON = Z_DON + Z_PON;
Z_OC = Z_DOC + Z_POC;

x = nan([size(P_C) nvars]);
for i = 1:nvars
    x(:,:,:,:,i) = eval(vars{i});
end

% If multiple trajectories selected then take averages
ntraj = length(traj);
if ntraj > 1
    x = squeeze(mean(x, ndims(x)-1, 'omitnan'));
end

% Integrate over depth
x = squeeze(sum(fixedParams.zwidth(:) .* permute(x, [2 1 3 4])));

x = x(:,:,strcmp(var, vars));

switch var
    case {'P_biovolume','Z_biovolume'}, x = 1e-9 .* x;
end

if logScale
   x = log10(x+1);
%    ensurePositive = false;
end

% if exist('minVal', 'var')
%     minVal = eval('minVal');
%     x(x < minVal) = nan;
% else
%     minVal = [];
% end

trophicLevel = var(1);

% Time-size grid for interpolation
[Size, Time] = ndgrid(abs(fixedParams.([trophicLevel 'Pdia'])), 1:fixedParams.nt);
[sizeGrid, timeGrid] = ndgrid(1:1:fixedParams.([trophicLevel 'Pdia_intervals'])(end), 1:fixedParams.nt);

if exist('xLim', 'var')
    sizeGrid = sizeGrid(:,xLim(1):xLim(2));
    timeGrid = timeGrid(:,xLim(1):xLim(2));
end
if exist('yLim', 'var')
    sizeGrid = sizeGrid(yLim(1):yLim(2),:);
    timeGrid = timeGrid(yLim(1):yLim(2),:);
end

% x = x .^ 0.1;
F = griddedInterpolant(Size, Time, x, smooth);
Fsmooth = F(sizeGrid, timeGrid);
if ensurePositive, Fsmooth(Fsmooth <= 0) = min(Fsmooth(Fsmooth > 0)); end
if exist('minVal', 'var')
    Fsmooth(Fsmooth < minVal) = nan;
end
% if ~isempty(minVal)
% %     Fsmooth(Fsmooth < minVal) = minVal;
%     Fsmooth(Fsmooth < minVal) = nan;
% end

if exist('levels', 'var')
    contourf(Fsmooth, levels)%, 'ShowText', 'on', 'LabelSpacing', 400)
else
   contourf(Fsmooth) 
end

set(gca, 'YScale', 'log')

cb = colorbar;

if logScale
    Ticks = get(cb, 'Ticks');
    tl = round(-1 + 10 .^ Ticks([1,end]), 1, 'significant');
    ltl = log10(tl+1);
    Ticks = linspace(ltl(1), ltl(2), length(Ticks));
    TickLabels = round(-1 + 10 .^ Ticks, 2, 'significant');
    TickLabels = arrayfun(@(z) num2str(z), TickLabels(:), 'UniformOutput', false);
    set(cb, {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
end



% if logScale
%     Ticks_ = get(cb, 'Ticks');
%     Ticks = round(10 .^ Ticks_, 2, 'significant');
%     iten = ismember(Ticks,10 .^ (-6:1:15));
%     if sum(iten) > 1
%         iten(end) = true;
%         Ticks = Ticks(iten);
%     end
%     TickLabels = arrayfun(@(z) num2str(z), -1 + Ticks(:), 'UniformOutput', false);
%     Ticks = log10(Ticks);    
%     set(cb, {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
% end


% if logScale
%     Ticks = get(cb, 'Ticks');
%     Ticks = round(10 .^ Ticks, 2, 'significant');
%     iten = ismember(Ticks,10 .^ (-6:1:10));
%     if sum(iten) > 1
%         iten(end) = true;
%         Ticks = Ticks(iten);
%     end
%     TickLabels = arrayfun(@(z) num2str(z), Ticks(:), 'UniformOutput', false);
%     Ticks = log10(Ticks);
%     set(cb, {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
% end
% 
% if exist('levels', 'var')
%     levels_ = levels(1,:);
%     levels_(1) = [];
%     set(cb, 'Ticks', levels_)
%     TickLabels = arrayfun(@(z) num2str(z), levels_(:), 'UniformOutput', false);
%     set(cb, 'TickLabels', TickLabels)
% end
% 

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



        