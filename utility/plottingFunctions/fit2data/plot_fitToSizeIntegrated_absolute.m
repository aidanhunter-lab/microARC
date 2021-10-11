function plot_fitToSizeIntegrated_absolute(xvar, trophicGroup, Data, modData, varargin)

% Compare nutrient (scalar) data to modelled equivalents -- data is grouped
% by depth and displayed as boxplots

%% Plotting options

extractVarargin(varargin)

if ~exist('colDat', 'var')
    colDat = [0, 0, 0];
end
if ~exist('colMod', 'var')
    colMod = [0, 1, 0];
end
if ~exist('logScale', 'var')
    logScale = 'semilogx'; % by default use log scale on x-axis
end
if ~exist('waterOrigin', 'var')
     % by default do not distinguish trajectories by water origin... this 
     % likely wont matter as data is probably already filtered by water origin
    waterOrigin = 'all';
end
if ~exist('connectors', 'var')
    connectors = false; % by default do not connect data to associated modelled points
end
if ~exist('meanOnly', 'var')
    meanOnly = false;
end
if ~exist('includeLegend', 'var')
    includeLegend = false;
end
if ~exist('legendLocation', 'var')
    legendLocation = 'northwest';
end


%% Dependencies on data type (xvar)

switch xvar
    case 'CellConc'
        Var = 'cell concentration';
        ylab = 'cell conc. (cells m^{-3})';
        ylab_rel = 'relative cell conc.';
        ylab_tot = 'total cell conc. (cells m^{-3})';
    case 'BioVol'
        Var = 'bio-volume';
        ylab = 'bio-volume (\mum^3 m^{-3})';
        ylab_rel = 'relative bio-volume';
        ylab_tot = 'total bio-volume (\mum^3 m^{-3})';
end


%% Create plot
ind = strcmp(Data.size.dataBinned.Variable, xvar); % index variable
ind = ind & strcmp(Data.size.dataBinned.trophicLevel, trophicGroup); % index trophic group
switch waterOrigin % index sample water origin
    case 'all'
%         ind = ind;
    case unique(Data.size.dataBinned.waterMass)
        ind = ind & strcmp(Data.size.dataBinned.waterMass, waterOrigin);
    otherwise
        warning('Optional argument "waterOrigin" does not correspond to values in the data. Not necessarily a problem, but check for mistakes...')
end

x = unique(Data.size.dataBinned.size(ind)); % cell diameter
ydat = Data.size.dataBinned.Value(ind);  % observed data
ymod = modData.size.Value(ind,:);        % modelled equivalent
% Average modelled values over sample events (& trajectories). We can plot
% values for all events, or just the means.
ym = mean(ymod(:,:), 2);
switch meanOnly
    case true
        ymod = ym;
        modMarkerFaceColor = colMod;
    case false
        modMarkerFaceColor = [1 1 1];
end

plotFun = @plot;
switch logScale
    case 'loglog'
        plotFun = @loglog;
    case 'semilogx'
        plotFun = @semilogx;
end

switch connectors, case true
    % Draw lines connecting data to mean modelled values
    for i = 1:length(ym)
        plotFun([x(i) x(i)], [ydat(i) ym(i)], 'Color', [0.85, 0.85, 0.85])
        if i == 1, hold on; end
    end
end

for i = 1:size(ymod, 2)
    pmod = plotFun(x, ymod(:,i), 'o', 'MarkerEdgeColor', colMod, ...
        'MarkerFaceColor', modMarkerFaceColor);
    if i == 1 && ~connectors, hold on; end
end
    
pdat = plotFun(x, ydat, 'o', 'Color', colDat);

hold off

xlabel('ESD (\mum)')
ylabel(ylab)

switch includeLegend, case true
    legend([pdat, pmod], 'data', 'model', 'Location', legendLocation)
end
