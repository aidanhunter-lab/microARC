function plot_fitToNutrient_event(xvar, Data, modData, varargin)

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
if ~exist('xLab', 'var')
    autoLabel = true;
end
if ~exist('logScale', 'var')
    logScale = false; % natural scale by default
end
if ~exist('standardised', 'var')
    standardised = true; % by default plot shows standardised data
end
if logScale && standardised
    logScale = false;
    warning('Incompatible options: "logScale" & "standardised" were both set as "true". Option "logScale" reset to "false".')
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
            xUnit = 'mg chl m^{-3}';
    end
    switch logScale, case true
        xUnit = ['log_{10}(' xUnit ')'];
    end
    switch standardised, case false
        xLab = [xLab ' (' xUnit ')'];
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
switch logScale, case true
    yobs = log10(yobs);
    ymod = log10(ymod);
end

event = Data.scalar.Event(ind);
uevent = unique(event);
nevent = length(uevent);

% Create matrix structure for boxplot -- stores all data and model points
for i = nevent:-1:1
    ei = event == uevent(i);
    x = yobs(ei);
    allValues(1:numel(x),2*nevent-2*i+1) = x(:);
    x = ymod(ei,:);
    allValues(1:numel(x),2*nevent-2*i+2) = x(:);
end

allValues(~isnan(allValues) & allValues == 0) = nan;

% Make boxplot
bp = boxplot2(allValues, 'orientation', 'horizontal', 'barwidth', 0.7);
% Asthetics
for i = 1:size(allValues, 2)
    bp.out(i).Marker = '.';
end
for i = 1:size(allValues, 2) / 2
    bp.out(2*i-1).MarkerFaceColor = colDat;
    bp.out(2*i-1).MarkerEdgeColor = colDat;
    bp.out(2*i).MarkerFaceColor = colMod;
    bp.out(2*i).MarkerEdgeColor = colMod;
    bp.box(2*i-1).Color = colDat;
    bp.box(2*i).Color = colMod;
    bp.med(2*i-1).Color = [0 0 0];
    bp.med(2*i).Color = [0 0 0];
end


gc = gca;
tt = reshape(1:size(allValues,2), [2, 0.5 * size(allValues,2)]);
gc.YTick = mean(tt);
gc.YTickLabel = num2str(flip(uevent));
gc.TickLength = 0.5 * gc.TickLength;
xl = gc.XLim; yl = gc.YLim;
tt = 0.5 * (tt(2,1:end-1) + tt(1,2:end));

% Lines separating boxes
hold on
for i = 1:length(tt)
    line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
end
hold off

xlabel(xLab)


ylabel('event')
