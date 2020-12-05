function Data = prepareFittingData(obsDir, varargin)

% Load and clean fitting data, output as struct

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Inorganic nutrient

nutrientObsFile = 'PS99_2_nutrients_rev.csv';

dat_nut = readtable(fullfile(obsDir,nutrientObsFile), 'Format', 'auto');
date = split(dat_nut.Date, 'T');
dat_nut.Date = date(:,1);
dat_nut.Time = date(:,2);
dat_nut = movevars(dat_nut, 'Time', 'After', 'Date');
% dat_nut.Year = year(dat_nut.Date); % this 'year' function apparently causes problems on some MatLab versions
dat_nut.Year = cellfun(@(x) str2double(x), ...
    cellfun(@(x) x(1:4), dat_nut.Date, 'UniformOutput', false));
dat_nut = movevars(dat_nut, 'Year', 'After', 'Time');
dat_nut.Yearday = yearday(datenum(dat_nut.Date));
dat_nut = movevars(dat_nut, 'Yearday', 'After', 'Year');
% filter out data columns
% dat_nut = dat_nut(:,{'Year','Yearday','Event','Event_2','Latitude', ... 
%     'Longitude','Depth','NO3_NO2_corrected','Flag_NO3_NO2_c'});
dat_nut = dat_nut(:,{'Year','Yearday','Event','Event_2','Latitude', ... 
    'Longitude','Depth','NO3_NO2_corrected','Flag_NO3_NO2_c', ... 
    'PO4_corrected','Flag_PO4_c','Si_OH4_corrected','Flag_Si_OH4_c'});

% remove rows flagged as (potentially) poor quality measurements
dat_nut.Properties.VariableNames(contains( ... 
    dat_nut.Properties.VariableNames, 'corrected')) = {'N','P','Si'};
dat_nut.Properties.VariableNames(contains( ... 
    dat_nut.Properties.VariableNames, 'Flag')) = {'Flag_N','Flag_P','Flag_Si'};
dat_nut.N(dat_nut.Flag_N ~= 0) = nan; dat_nut.Flag_N = [];
dat_nut.P(dat_nut.Flag_P ~= 0) = nan; dat_nut.Flag_P = [];
dat_nut.Si(dat_nut.Flag_Si ~= 0) = nan; dat_nut.Flag_Si = [];

% convert from short to long format
dat_nut = stack(dat_nut, {'N','P','Si'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat_nut.Variable = cellstr(dat_nut.Variable);
dat_nut(isnan(dat_nut.Value),:) = [];
dat_nut.Type = repmat({'inorganic'}, [height(dat_nut) 1]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% POM and Chl_a

OMObsFile = 'copy_data_Engel_etal2019.csv';

dat_OM = readtable(fullfile(obsDir,OMObsFile), 'Format', 'auto');
% The lat-long column entries are not formatted consistently across rows,
% and also contain degree symbols and measurements both in decimal-degrees 
% and degree-minutes. This needs fixed...
ind = ~cellfun('isempty', strfind(dat_OM.Latitude, char(176)));
degr = split(dat_OM.Latitude(ind),char(176)); % convert degree-minutes/seconds into degree-decimal
minute = split(degr(:,2),'.');
degr = cellfun(@(x)str2double(x),degr(:,1));
second = cellfun(@(x)str2double(x),minute(:,2));
minute = cellfun(@(x)str2double(x),minute(:,1));
lat = cellstr(num2str(degr + minute ./ 60 + second ./ 3600));
dat_OM.Latitude(ind) = lat;
dat_OM.Latitude = cellfun(@(x)str2double(x), dat_OM.Latitude);

dat_OM.Longitude = strrep(dat_OM.Longitude,'E00','');
dat_OM.Longitude = strrep(dat_OM.Longitude,'E0','');
dat_OM.Longitude = strrep(dat_OM.Longitude,'E','');
dat_OM.Longitude = strrep(dat_OM.Longitude,'W00','-');
dat_OM.Longitude = strrep(dat_OM.Longitude,'W0','-');
dat_OM.Longitude = strrep(dat_OM.Longitude,'W','-');
ind = ~cellfun('isempty', strfind(dat_OM.Longitude, char(176)));
degr = split(dat_OM.Longitude(ind),char(176));
minute = split(degr(:,2),'.');
degr = cellfun(@(x)str2double(x),degr(:,1));
second = cellfun(@(x)str2double(x),minute(:,2));
minute = cellfun(@(x)str2double(x),minute(:,1));
lon = cellstr(num2str(degr + minute ./ 60 + second ./ 3600));
dat_OM.Longitude(ind) = lon;
dat_OM.Longitude = cellfun(@(x)str2double(x), dat_OM.Longitude);
% Use only 2017 from this data set, because we have particle trajectories
% for this year.
dat_OM = dat_OM(dat_OM.Year == 2017,:);
dat_OM.Yearday = yearday(datenum(dat_OM.Date));
% filter out data columns
dat_OM = dat_OM(:,{'Year','Yearday','ARK_oder_PS','HG_stations','Station','EventNumber','Latitude','Longitude','Depth','chl_a','POC','PON'});
% Convert POM units from (mu g / L) to (m mol / m^3)
dat_OM.PON = dat_OM.PON ./ 14;
dat_OM.POC = dat_OM.POC ./ 12;
% % Convert chl units from mu g / L to mu g / m^3
% dat_OM.chl_a = dat_OM.chl_a * 1000;
% dat_OM.Properties.VariableNames{'chl_a'} = 'Chl';
% convert from short to long format
dat_OM = stack(dat_OM, {'chl_a','POC','PON'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat_OM.Variable = cellstr(dat_OM.Variable);
dat_OM(isnan(dat_OM.Value),:) = [];
dat_OM.Type = repmat({'organic'}, [height(dat_OM) 1]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Size spectra

% load data
sizeSpectraObsFile = 'S1-S4_spectra_Nbiomass.csv';
dat_size = readtable(fullfile(obsDir, sizeSpectraObsFile), 'Format', 'auto');

% A selection of cell size -> N conversions are contained in these data.
% Although the different conversions produce nitrogen size-spectra of
% similar shape, there is some substantial variation in magnitudes
% (especially at small and large cell sizes). So maybe we should
% use the measured cell densities in our cost function rather than the
% nitrogen densities...
varNames = dat_size.Properties.VariableNames;
varNames = varNames(contains(varNames, 'NperCell'));
ne = length(varNames);
NsizeEqs = cell(1,ne); % Store names of size -> nitrogen conversions
for i = 1:ne
    varSplit = strsplit(varNames{i}, '_');
    NsizeEqs{i} = varSplit{2};
end

% ESDmin = 1; % min/max sizes to retain from the data
% ESDmax = 200;
% if ESDmin < min(dat_size.ESD), ESDmin = min(dat_size.ESD); end
% if ESDmax > max(dat_size.ESD), ESDmax = max(dat_size.ESD); end

% Convert units of original data into units used in model (keep all
% elemental concentrations in mmol)
dat_size.cellDensity = 1e3 * dat_size.cellDensity; % cells/L -> cells/m^3
dat_size.cellDensitySD = 1e3 * dat_size.cellDensitySD; % cells/L -> cells/m^3
dat_size.cellDensitySE = 1e3 * dat_size.cellDensitySE; % cells/L -> cells/m^3
dat_size.cellDensityCImin = 1e3 * dat_size.cellDensityCImin; % cells/L -> cells/m^3
dat_size.cellDensityCImax = 1e3 * dat_size.cellDensityCImax; % cells/L -> cells/m^3
for i = 1:ne
    % nitrogen per cell: pg N / cell -> mmol N / cell
    dat_size.(['NperCell_' NsizeEqs{i}]) = ... 
        (1/14) * 1e-9 * dat_size.(['NperCell_' NsizeEqs{i}]);
end
for i = 1:ne
    % nitrogen density: mug N / L -> mmol N / m^3
    dat_size.(['Ndensity_' NsizeEqs{i}]) = ... 
        (1/14) * dat_size.(['Ndensity_' NsizeEqs{i}]);
end

% Take average values over the various size -> N conversions
dat_size.NperCell = mean(dat_size{:,strcat('NperCell_', NsizeEqs)}, 2);
dat_size.Ndensity = mean(dat_size{:,strcat('Ndensity_', NsizeEqs)}, 2);

% Remove unnecessary variables and rows from data
scenarios = unique(dat_size.scenario)'; % data collected from different cruises/years
N = length(scenarios);
TrophicLevel = 'autotroph'; % At present we're only modelling size spectra of autotrophs
% keepRows = strcmp(dat_size.trophicLevel, TrophicLevel) & ... 
%     dat_size.ESD >= ESDmin & dat_size.ESD <= ESDmax;
keepRows = strcmp(dat_size.trophicLevel, TrophicLevel);
dat_size = dat_size(keepRows,:);
keepVars = [{'scenario', 'season', 'regime', 'ESD', 'cellVolume', ...
    'cellDensity', 'cellDensitySD'}, strcat('Ndensity_', NsizeEqs), 'Ndensity'];
dat_size = dat_size(:,keepVars);

cols = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]]; % plotting colours for different scenarios

v = reshape(varargin, [2 0.5*length(varargin)]);

%~~~
% optional plots
plotNconcSpectra = ~isempty(v) && any(contains(v(1,:),'plotNconcSpectra')) && v{2,strcmp(v(1,:), 'plotNconcSpectra')};
plotCellConcSpectra = ~isempty(v) && any(contains(v(1,:),'plotCellConcSpectra')) && v{2,strcmp(v(1,:), 'plotCellConcSpectra')};
if plotNconcSpectra || plotCellConcSpectra
    figure
    if plotNconcSpectra && plotCellConcSpectra
        subplot(1,2,1)
        ii = strcmp(dat_size.scenario, scenarios(1));
        loglog(dat_size.ESD(ii), dat_size.Ndensity(ii), 'Color', cols(1,:))
        hold on
        for i = 2:N
            ii = strcmp(dat_size.scenario, scenarios(i));
            loglog(dat_size.ESD(ii), dat_size.Ndensity(ii), 'Color', cols(i,:))
        end
        gc = gca;
        gc.YLim(1) = 1e-6;
        gc.YLim(2) = 2 * max(dat_size.Ndensity);
        xlabel('ESD (\mum)')
        ylabel('Nitrogen density (mmol N m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$)', 'Interpreter', 'latex')
        %         title('Size spectra data')
        for i = N:-1:1
            j = N-i+1;
            base = 5; space = 5;
            yt = space^(j-1) * base * gc.YLim(1);
            text(gc.XLim(1) * (1 + 0.3 * 10), yt, scenarios{i});
            line([gc.XLim(1) * (1 + 0.1 * 10) gc.XLim(1) * (1 + 0.2 * 10)], [yt yt], 'Color', cols(i,:))
        end
        hold off
        subplot(1,2,2)
        ii = strcmp(dat_size.scenario, scenarios(1));
        loglog(dat_size.ESD(ii), dat_size.cellDensity(ii), 'Color', cols(1,:))
        hold on
        for i = 2:N
            ii = strcmp(dat_size.scenario, scenarios(i));
            loglog(dat_size.ESD(ii), dat_size.cellDensity(ii), 'Color', cols(i,:))
        end
        gc = gca;
        gc.YLim(1) = 1e1;
        gc.YLim(2) = 2 * max(dat_size.cellDensity);
        xlabel('ESD (\mum)')
        ylabel('cell conc. density (cells m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$)', 'Interpreter', 'latex')
        %         title('Size spectra data')
        for i = N:-1:1
            j = N-i+1;
            base = 5; space = 5;
            yt = space^(j-1) * base * gc.YLim(1);
            text(gc.XLim(1) * (1 + 0.3 * 10), yt, scenarios{i});
            line([gc.XLim(1) * (1 + 0.1 * 10) gc.XLim(1) * (1 + 0.2 * 10)], [yt yt], 'Color', cols(i,:))
        end
        hold off
        sgtitle('Size spectra data')
    else        
        if plotNconcSpectra
            ii = strcmp(dat_size.scenario, scenarios(1));
            loglog(dat_size.ESD(ii), dat_size.Ndensity(ii), 'Color', cols(1,:))
            hold on
            for i = 2:N
                ii = strcmp(dat_size.scenario, scenarios(i));
                loglog(dat_size.ESD(ii), dat_size.Ndensity(ii), 'Color', cols(i,:))
            end
            gc = gca;
            gc.YLim(1) = 1e-6;
            gc.YLim(2) = 2 * max(dat_size.Ndensity);
            xlabel('ESD (\mum)')
            ylabel('Nitrogen density (mmol N m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$)', 'Interpreter', 'latex')
            title('Size spectra data')
            for i = N:-1:1
                j = N-i+1;
                base = 5; space = 5;
                yt = space^(j-1) * base * gc.YLim(1);
                text(gc.XLim(1) * (1 + 0.3 * 10), yt, scenarios{i});
                line([gc.XLim(1) * (1 + 0.1 * 10) gc.XLim(1) * (1 + 0.2 * 10)], [yt yt], 'Color', cols(i,:))
            end
            hold off
        end        
        if plotCellConcSpectra
            ii = strcmp(dat_size.scenario, scenarios(1));
            loglog(dat_size.ESD(ii), dat_size.cellDensity(ii), 'Color', cols(1,:))
            hold on
            for i = 2:N
                ii = strcmp(dat_size.scenario, scenarios(i));
                loglog(dat_size.ESD(ii), dat_size.cellDensity(ii), 'Color', cols(i,:))
            end
            gc = gca;
            gc.YLim(1) = 1e1;
            gc.YLim(2) = 2 * max(dat_size.cellDensity);
            xlabel('ESD (\mum)')
            ylabel('cell conc. density (cells m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$)', 'Interpreter', 'latex')
            title('Size spectra data')
            for i = N:-1:1
                j = N-i+1;
                base = 5; space = 5;
                yt = space^(j-1) * base * gc.YLim(1);
                text(gc.XLim(1) * (1 + 0.3 * 10), yt, scenarios{i});
                line([gc.XLim(1) * (1 + 0.1 * 10) gc.XLim(1) * (1 + 0.2 * 10)], [yt yt], 'Color', cols(i,:))
            end
            hold off
        end
    end
end
%~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Merge data sets: group all the 'scalar' data into single table and store
% the size data separately.

% Match column names
dat_nut.Properties.VariableNames{'Event_2'} = 'Station';
dat_nut.Properties.VariableNames{'Event'} = 'EventLabel';
dat_nut = movevars(dat_nut, 'EventLabel', 'After', 'Station');
dat_OM = removevars(dat_OM, {'ARK_oder_PS','Station'});
dat_OM.Properties.VariableNames{'EventNumber'} = 'EventLabel';
dat_OM.Properties.VariableNames{'HG_stations'} = 'Station';

% Station names used in each data set are written inconsistently
% disp(unique(dat_nut.Station)')
% disp(unique(dat_OM.Station)')
dat_OM.Station = strrep(dat_OM.Station,' shallow','');
dat_OM.Station = strrep(dat_OM.Station,'I ','I');
dat_OM.Station = strrep(dat_OM.Station,' ','_');
dat_nut.Station = strrep(dat_nut.Station,'S_3','S3');

% Event format is inconsistent
dat_nut.EventLabel = strrep(dat_nut.EventLabel, '/', '_');
dat_nut.EventLabel = strrep(dat_nut.EventLabel, '-', '_');
dat_OM.EventLabel = strrep(dat_OM.EventLabel, '/', '_');
dat_OM.EventLabel = strrep(dat_OM.EventLabel, '-', '_');

% Match dates & depths to scenarios in size data
scenarioInfo.S1.dateFirst = datetime(2016,6,23);
scenarioInfo.S1.dateLast = datetime(2016,7,16);
scenarioInfo.S2.dateFirst = datetime(2017,7,23);
scenarioInfo.S2.dateLast = datetime(2017,8,19);
scenarioInfo.S3.dateFirst = datetime(2018,7,10);
scenarioInfo.S3.dateLast = datetime(2018,8,3);
scenarioInfo.S4.dateFirst = datetime(2018,9,16);
scenarioInfo.S4.dateLast = datetime(2018,10,12);
scenarioInfo.S1.depth = [10 30];
scenarioInfo.S2.depth = [10 40];
scenarioInfo.S3.depth = [5 43];
scenarioInfo.S4.depth = [5 34];

dat_size.Year = nan(height(dat_size), 1);
% dat_size.Yearday = nan(height(dat_size), 1);
dat_size.DepthMin = nan(height(dat_size), 1);
dat_size.DepthMax = nan(height(dat_size), 1);
fields = fieldnames(scenarioInfo);
for i = 1:length(fields)
    ind = strcmp(dat_size.scenario,fields{i});
    indn = sum(ind);
    [yr, ~] = datevec(scenarioInfo.(fields{i}).dateFirst);
    yrDayFirst = yearday(datenum(scenarioInfo.(fields{i}).dateFirst));
    yrDayLast = yearday(datenum(scenarioInfo.(fields{i}).dateLast));
    dat_size.Year(ind) = repmat(yr, [indn 1]);
    dat_size.YeardayFirst(ind) = repmat(yrDayFirst, [indn 1]);    
    dat_size.YeardayLast(ind) = repmat(yrDayLast, [indn 1]);    
    dat_size.DepthMin(ind) = scenarioInfo.(fields{i}).depth(1);
    dat_size.DepthMax(ind) = scenarioInfo.(fields{i}).depth(2);    
end

dat_size = movevars(dat_size, {'Year', 'YeardayFirst', 'YeardayLast'}, 'Before', 'season');
dat_size = movevars(dat_size, {'DepthMin', 'DepthMax'}, 'After', 'regime');


% Merge nutrient data (inorganic & organic)
dat = [dat_nut; dat_OM];

% Remove NaNs & rows with non-positive measures
dat(dat.Value <= 0 | isnan(dat.Value),:) = [];

% Index each unique sampling event
EventLabel = unique(dat.EventLabel, 'stable');
Event = (1:length(EventLabel))';
events = table(Event, EventLabel);
dat = join(dat, events);

% % As nitrogen is the only modelled inorganic nutrient (this may change in 
% % future), omit measures of other nutrients.
% remove = strcmp(dat.Type, 'inorganic') & ~strcmp(dat.Variable, 'N');
% dat(remove,:) = [];

% Store fitting data as a struct
scalarData = table2struct(dat, 'ToScalar', true);
scalarData.Variable = cellstr(scalarData.Variable);
scalarData.nSamples = size(scalarData.Value,1);
scalarData.nEvents = length(unique(scalarData.Event));

sizeData = table2struct(dat_size, 'ToScalar', true);
% sizeData.nSamples = size(sizeData.Year,1);

% Include fields specifying which of these data are used in the cost
% function
% scalarData.obsInCostFunction = {'N', 'chl_a', 'PON', 'POC'};
% scalarData.inCostFunction = ismember(scalarData.Variable, ... 
%     scalarData.obsInCostFunction);

% sizeData.obsInCostFunction = {'NConc'}; % maybe move this to function that selects size class intervals...
% sizeData.dataBinned.inCostFunction = ismember(sizeData.dataBinned.Variable, ... 
%     sizeData.obsInCostFunction);

% % Omit any events where data used in the cost function was not collected
% fields = fieldnames(scalarData);
% uev = unique(scalarData.Event);
% for i = 1:length(uev)
%     ind = scalarData.Event == uev(i);
%     x = unique(scalarData.Variable(ind));
%     y = ismember(x, scalarData.obsInCostFunction);
%     if ~any(y)
%         for j = 1:length(fields)
%             if size(scalarData.(fields{j}), 1) > 1
%                 scalarData.(fields{j}) = scalarData.(fields{j})(~ind);
%             end
%         end
%     end
% end
% if length(unique(scalarData.Event)) < scalarData.nEvents
%     scalarData.nSamples = size(scalarData.Value,1);    
%     scalarData.nEvents = length(unique(scalarData.Event));    
%     oldEventLabs = unique(scalarData.Event);
%     newEventLabs = (1:length(oldEventLabs))';    
%     x = table(oldEventLabs, newEventLabs);
%     for i = 1:height(x)
%         scalarData.Event(scalarData.Event == x.oldEventLabs(i)) = x.newEventLabs(i);
%     end
% end



Data.scalar = scalarData;
Data.size = sizeData;



end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~
% Functions
%~~~~~~~~~~

% Use stationary and inflection points of size spectra to generate size
% class intervals

function y = stationaryPoints(x)
% Return positions of turning points in x (and also endpoints).
% Only accepts vector arguments... more work needed for matrix inputs
% because the output will need to be a cell array where each component
% holds a vector of arbitrary length...
xdim0 = size(x);
vec = isvector(x);
if vec
    x_ = x(:);
    n = length(x);
else
    x_ = x;
end
xdim = size(x_);
y = diff([x_(1,:); x_]);
y = y > 0;
y = (y(1:end-1,:) & ~y(2:end,:)) | (~y(1:end-1,:) & y(2:end,:));
y = [y; zeros(1, xdim(2))];
[r, c] = find(y);
if vec
    y = r;
    y = [1; y; n];
else
    yy = cell(1,xdim(2));
    for i = 1:xdim(2)
        yy{i} = r(c == i);
        yy{i} = [1; yy{i}; dim(1)];
    end
    y = yy;
end
if all(xdim == flip(xdim0)), y = y'; end
end

function y = inflectionPoints(x)
% Return positions of inflection points in x (and also endpoints)
xdim0 = size(x);
vec = isvector(x);
if vec
    x_ = x(:);
    n = length(x);
else
    x_ = x;
end
xdim = size(x_);
y = diff(diff(x_));
y = [y(1,:); y; y(end,:)];
y = y > 0;
y = (y(1:end-1,:) & ~y(2:end,:)) | (~y(1:end-1,:) & y(2:end,:));
y = [y; zeros(1, xdim(2))];
[r, c] = find(y);
if vec
    y = r;
    y = [1; y; n];
else
    yy = cell(1,xdim(2));
    for i = 1:xdim(2)
        yy{i} = r(c == i);
        yy{i} = [1; yy{i}; xdim(1)];
    end
    y = yy;
end
if all(xdim == flip(xdim0)), y = y'; end
end


function y = partitionSizeSpectra2(x, varargin)

% Input size spectra vectors, x.
% Choose data features in varargin.
% Returns the positions, y, in x where those features occur, which can be
% used to partition the size spectra vectors into size class bins.
if isempty(varargin)
    % By default use stationary points
    Type = 'stationary';
    % By default use stationary and inflecion points, not the midpoints
    % between them
    useMidpoints = false;
else
    % Currently have options to use stationary and/or inflection points
    i = strcmp(varargin, 'stationaryPoints');
    if any(i)
        if varargin{find(i)+1}
            Type = 'stationary';
        end
    end
    i = strcmp(varargin, 'inflectionPoints');
    if any(i)
        if varargin{find(i)+1}
            if ~exist('Type', 'var')
                Type = 'inflection';
            else
                Type = 'both';
            end
        end
    end
    i = strcmp(varargin, 'useMidpoints');
    if any(i)
        useMidpoints = varargin{find(i)+1};
    end
end
if ~exist('Type', 'var')
    y = [];
    warning('Check validity of optional arguments in partitionSizeSpectra.m')
else
    xdim = size(x);
    if length(xdim) ~= 2 || ~any(xdim == 1)
        warning('partitionSizeSpectra only accepts vector arguments')
    else
        switch Type
            case 'stationary'
                y = stationaryPoints(x);
                if useMidpoints
                    y = round([y(1,:); 0.5 * (y(1:end-1,:) + y(2:end,:)); y(end,:)]);
                end                
            case 'inflection'
                y = inflectionPoints(x);
                if useMidpoints
                    y = round([y(1,:); 0.5 * (y(1:end-1,:) + y(2:end,:)); y(end,:)]);
                end
            case 'both'
                y1 = stationaryPoints(x);
                y2 = inflectionPoints(x);
                if useMidpoints
                    y1 = round([y1(1,:); 0.5 * (y1(1:end-1,:) + y1(2:end,:)); y1(end,:)]);
                    y2 = round([y2(1,:); 0.5 * (y2(1:end-1,:) + y2(2:end,:)); y2(end,:)]);
                end
                y = unique([y1(:); y2(:)]);
        end
    end
end
end






% function y = partitionSizeSpectra2(x, varargin)
% 
% % Input size spectra vectors, x.
% % Choose data features in varargin.
% % Returns the positions, y, in x where those features occur, which can be
% % used to partition the size spectra vectors into size class bins.
% 
% if isempty(varargin)
%     % By default use stationary points
%     Type = 'stationary';
% else
%     % Currently have options to use stationary and/or inflection points
%     i = strcmp(varargin, 'stationaryPoints');
%     if any(i)
%         if varargin{find(i)+1}
%             Type = 'stationary';
%         end
%     end
%     i = strcmp(varargin, 'inflectionPoints');
%     if any(i)
%         if varargin{find(i)+1}
%             if ~exist('Type', 'var')
%                 Type = 'inflection';
%             else
%                 Type = 'both';
%             end
%         end
%     end
% end
% 
% if ~exist('Type', 'var')
%     y = [];
%     warning('Check validity of optional arguments in partitionSizeSpectra.m')
% else
%     
%     n = length(x);
%     xdim = size(x);
%     
%     if length(xdim) ~= 2 || ~any(xdim == 1)
%         warning('partitionSizeSpectra only accepts vector arguments')
%     else
%         switch Type
%             case 'stationary'
%                 y = diff(x(:));
%                 y = [y(1); y];
%                 y = y > 0;
%                 y = (y(1:n-1) & ~y(2:n)) | (~y(1:n-1) & y(2:n));
%                 y = [0; y];
%                 y = [1; find(y); n]; % endpoints and turning points
%                 if xdim(1) == 1
%                     y = y';
%                 end
%                 
%             case 'inflection'
%                 y = diff(diff(x(:)));
%                 y = [y(1); y(1); y];
%                 y = y > 0;
%                 y = (y(1:n-1) & ~y(2:n)) | (~y(1:n-1) & y(2:n));
%                 y = [0; y];
%                 y = [1; find(y); n]; % endpoints and inflection points
%                 if xdim(1) == 1
%                     y = y';
%                 end
%                 
%             case 'both'
%                 yd = diff(x(:));
%                 y = [yd(1); yd];
%                 y = y > 0;
%                 y = (y(1:n-1) & ~y(2:n)) | (~y(1:n-1) & y(2:n));
%                 y = [0; y];
%                 y0 = [1; find(y); n]; % endpoints and turning points
%                 
%                 yd = diff(yd);
%                 y = [yd(1); yd(1); yd];
%                 y = y > 0;
%                 y = (y(1:n-1) & ~y(2:n)) | (~y(1:n-1) & y(2:n));
%                 y = [0; y];
%                 y = [1; find(y); n]; % endpoints and inflection points
%                 y = unique([y; y0]); % endpoints, turning points, and inflection points
%         end
%     end
% end
