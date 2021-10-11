function plot_fitToNutrient_event2(xvar, Data, modData, varargin)

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
if ~exist('omitSingles', 'var')
    omitSingles = false;
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
switch omitSingles
    case false
        for i = nevent:-1:1
            ei = event == uevent(i);
            x = yobs(ei);
            allValues(1:numel(x),2*nevent-2*i+1) = x(:);
            x = ymod(ei,:);
            allValues(1:numel(x),2*nevent-2*i+2) = x(:);
        end
    case true
        omit_i = false(nevent, 1);
        omit_ind = false(size(yobs));
        for i = nevent:-1:1
            ei = event == uevent(i);
            omitSingle = sum(ei) == 1;
            if omitSingle
                omit_i(nevent-i+1) = true;
                omit_ind = omit_ind | ei;
            end
            x = yobs(ei);
            allValues(1:numel(x),2*nevent-2*i+1) = x(:);
            x = ymod(ei,:);
            allValues(1:numel(x),2*nevent-2*i+2) = x(:);
        end
        allValues(:,repelem(omit_i, 2)) = [];
        event(omit_ind) = [];
        uevent = unique(event);
        nevent = length(uevent);
        
    case 'merge'
        
        merge_i = false(nevent, 1);
        merge_ind = false(size(yobs));
        for i = nevent:-1:1
            ei = event == uevent(i);
            mergeSingle = sum(ei) == 1;
            if mergeSingle
                merge_i(nevent-i+1) = true;
                merge_ind = merge_ind | ei;
            end
            x = yobs(ei);
            allValues(1:numel(x),2*nevent-2*i+1) = x(:);
            x = ymod(ei,:);
            allValues(1:numel(x),2*nevent-2*i+2) = x(:);
        end
        mergeEvents = uevent(flip(merge_i));
        plotEvents = uevent(flip(~merge_i));
        if ~isempty(mergeEvents)
            allValues = [allValues; zeros(1, size(allValues, 2))];
        end
        for i = 1:length(mergeEvents)
            md = mergeEvents(i);
            pd = plotEvents(abs(plotEvents - md) == min(abs(plotEvents - md)));
            pd = max(pd);
            allValues(end,repelem(flip(uevent) == pd, 2)) = ...
                allValues(1,repelem(flip(uevent) == md, 2));
        end
        allValues(:,repelem(merge_i, 2)) = [];
        event(merge_ind) = [];
        uevent = unique(event);
        nevent = length(uevent);
end


allValues(~isnan(allValues) & allValues == 0) = nan;

% Make boxplot

xx = reshape(1:3 * size(allValues, 2) / 2, 3, []);
xx = xx(1:2,:);

Cols = repmat([colDat; colMod], nevent, 1);

bwLo = 0.15;
bwHi = 0.9;
bw = sum(~isnan(allValues)) .^ 0.25;
bw = (bw - min(bw)) ./ (max(bw) - min(bw));
bw = bwLo + (bwHi-bwLo) * bw;
MarkerStyle = '.';

for i = 1:size(allValues, 2)
    boxchart(xx(i)*ones(size(allValues(:,i))), allValues(:,i), ...
        'orientation', 'horizontal', 'MarkerStyle', MarkerStyle, ...
        'MarkerColor', Cols(i,:), 'BoxFaceColor', Cols(i,:), 'BoxWidth', bw(i))
    if i == 1, hold on; end
end

set(gca, {'YTick', 'YTickLabel'}, {mean(xx), num2str(flip(uevent))})
set(gca, 'TickLength', 0.5 * get(gca, 'TickLength'))


% Lines separating boxes
hold on
tt = 0.5 * (xx(1,2:end) + xx(end,1:end-1));
xl = get(gca, 'XLim');
for i = 1:length(tt)
    line(xl, [tt(i) tt(i)], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
end
hold off

xlabel(xLab)

ylabel('event')

set(gca, 'XLim', xl)
set(gca, 'YLim', [0 max(xx(:))+1])


