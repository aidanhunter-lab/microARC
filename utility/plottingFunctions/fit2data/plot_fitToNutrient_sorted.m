function plot_fitToNutrient_sorted(xvar, Data, modData, varargin)

% Compare nutrient (scalar) data to modelled equivalents -- data is not
% grouped by any covariate, and is sorted into ascending order.

%% Plotting options

extractVarargin(varargin)

if ~exist('colDat', 'var')
    colDat = [0, 0, 0];
end
if ~exist('colMod', 'var')
    colMod = [0, 1, 0];
end
if ~exist('xLab', 'var')
    autoLabel = true;
end
if ~exist('standardised', 'var')
    standardised = true; % by default plot shows standardised data
end
if ~exist('connectors', 'var')
    connectors = false; % by default do not connect data to associated modelled points
end


%% Dependencies on data type (xvar)

if autoLabel
    switch xvar
        % Plot labels
        case 'N'
            xLab = 'DIN';
            xUnit = 'mmol N m^{-3}';
        case 'PON'
            xLab = 'PON';
            xUnit = 'mmol N m^{-3}';
        case 'POC'
            xLab = 'POC';
            xUnit = 'mmol C m^{-3}';
        case 'chl_a'
            xLab = 'chl a';
            xUnit = 'mg chl a m^{-3}';
    end
end

%% Create plot

ind = strcmp(Data.scalar.Variable, xvar);
switch standardised
    case true
        yobs = Data.scalar.scaled_Value(ind);
        ymod = modData.scalar.scaled_Value(ind,:);
    case false
        yobs = Data.scalar.Value(ind);
        ymod = modData.scalar.Value(ind,:);
end

n = length(yobs);
p = (1:n) ./ n;

[yobs_sort, o] = sort(yobs);
ymod_sort = ymod(o,:);

switch connectors, case true
    for ij = 1:n
        % tracer-lines connecting data and associated modelled points
        plot([yobs_sort(ij,1), ymod_sort(ij,1)], [p(ij), p(ij)], 'Color', [0.85, 0.85, 0.85])
        if ij == 1, hold on; end
    end
end

plot(yobs_sort, p, 'Color', colDat, 'LineWidth', 2) % plot data as line

switch connectors, case false, hold on; end

scatter(ymod_sort, p, 'MarkerEdgeColor', colMod) % plot modelled equivalents as points

hold off

switch standardised
    case true
        xlabel(xLab)
    case false
        xlabel([xLab ' (' xUnit ')'])
end

gc = gca;
gc.YTickLabel = [];
