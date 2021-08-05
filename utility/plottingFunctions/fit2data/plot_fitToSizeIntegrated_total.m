function plot_fitToSizeIntegrated_total(xvar, trophicGroup, Data, modData, varargin)

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
    logScale = 'none'; % by default no transformed scale
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

ydat = Data.size.dataBinned.Value(ind);  % observed data
ymod = modData.size.Value(ind,:);        % modelled equivalent
ydat = sum(ydat);
ymod = sum(ymod);


% Average modelled values over sample events (& trajectories). We can plot
% values for all events, or just the means.
ym = mean(ymod, 2);
switch meanOnly
    case true
        ymod = ym;
        modMarkerFaceColor = colMod;
    case false
        modMarkerFaceColor = [1 1 1];
end

switch connectors, case true
    % Draw lines connecting data to mean modelled values
    for i = 1:length(ym)
        plot([0.5 0.5], [ydat(i) ym(i)], 'Color', [0.85, 0.85, 0.85])
        if i == 1, hold on; end
    end
end

pmod = plot(repmat(0.5, [1 length(ymod)]), ymod, 'o', 'MarkerFaceColor', modMarkerFaceColor, ...
    'MarkerEdgeColor', colMod);
if ~connectors, hold on; end

pdat = plot(0.5, ydat, 'o', 'MarkerEdgeColor', colDat);

gc = gca;
yl = gc.YLim;
yl(1) = 0;
yl(2) = 1.5 .* yl(2);
gc.YLim = yl;
gc.XLim = [0 1];
gc.XTick = [];

hold off

ylabel(ylab_tot)

switch includeLegend, case true
    legend([pdat, pmod], 'data', 'model', 'Location', legendLocation)
end
