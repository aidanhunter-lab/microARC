function Data = prepareFittingData(obsDir, varargin)

% Include time-of-day in data output

% Load and clean fitting data, output as struct


%% Inorganic nutrient

%~~~~~~~~~~~~
% PS99 cruise
%~~~~~~~~~~~~

nutrientObsFile = 'PS99_2_nutrients_rev.csv';

dat = readtable(fullfile(obsDir,nutrientObsFile), 'Format', 'auto');
date = split(dat.Date, 'T');
dat.Date = date(:,1);
dat.Time = date(:,2);
dat.Time = cellfun(@(x) x(1:5), dat.Time, 'UniformOutput', false); % omit the seconds
dat.t = datenum(dat.Date);
[dat.Year, ~] = datevec(dat.t);
dat.Yearday = yearday(dat.t);
% filter out columns
dat = dat(:,{'Event','t','Date','Time','Year','Yearday','Latitude', ... 
    'Longitude','Depth','NO3_NO2_corrected','Flag_NO3_NO2_c', ... 
    'PO4_corrected','Flag_PO4_c','Si_OH4_corrected','Flag_Si_OH4_c'});
% remove rows flagged as (potentially) poor quality measurements
dat.Properties.VariableNames(contains( ... 
    dat.Properties.VariableNames, 'corrected')) = {'N','P','Si'};
dat.Properties.VariableNames(contains( ... 
    dat.Properties.VariableNames, 'Flag')) = {'Flag_N','Flag_P','Flag_Si'};
dat.N(dat.Flag_N ~= 0) = nan; dat.Flag_N = [];
dat.P(dat.Flag_P ~= 0) = nan; dat.Flag_P = [];
dat.Si(dat.Flag_Si ~= 0) = nan; dat.Flag_Si = [];
% convert from short to long format
dat = stack(dat, {'N','P','Si'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat.Variable = cellstr(dat.Variable);
dat(isnan(dat.Value),:) = [];
dat.Type = repmat({'inorganic'}, [height(dat) 1]);
% reformat event label to include leading zero on sample number
dat.Properties.VariableNames{'Event'} = 'Label';
x = split(dat.Label, '-');
x2 = x(:,2);
for i = 1:length(x2)
    if length(x2{i}) == 1
        x2{i} = ['0' x2{i}];
    end
end
dat.Label = join([x(:,1), x2], '-');

dat_nut = dat; % nutrient data suffixed with '_nut'


%~~~~~~~~~~~~~
% PS114 cruise
%~~~~~~~~~~~~~

nutrientObsFile = 'PS114_Nuts_Final_PANGAEA_v1.csv';
stationList = 'PS114_station_list.csv';

dat = readtable(fullfile(obsDir,nutrientObsFile), 'Format', 'auto');
stations = readtable(fullfile(obsDir,stationList), 'Format', 'auto');
% There is no date provided in dat(!), so need to combine table of
% sampling station data with table of measurements... but the formating is
% inconsistent...
dat.Label = join([dat.Cruise, dat.Station], '/');
dat.Label = strrep(dat.Label, '_', '-');
stations.Properties.VariableNames{'Station'} = 'Label';
% remove underway measures
stations = stations(~contains(stations.Label, 'Underway'),:);
% reformat
stations.Label = strrep(stations.Label, '_', '/');
% include leading zeros
x = split(stations.Label, '/');
x2 = x(:,2);
x_ = split(x2, '-');
x_1 = x_(:,1);
x_2 = x_(:,2);
for i = 1:length(x_1)
    n1 = length(x_1{i});
    n2 = length(x_2{i});
    if n1 == 1, x_1{i} = ['00' x_1{i}]; end
    if n1 == 2, x_1{i} = ['0' x_1{i}]; end
    if n2 == 1, x_2{i} = ['0' x_2{i}]; end
end
x_ = join([x_1, x_2], '-');
stations.Label = join([x(:,1), x_], '/');

% extract event dates/times
stations = stations(:,{'Label', 'Date', 'Time', 'Action'});
stations = stations(contains(stations.Action, 'start'),:); % use event starting times
stations.Action = [];
% some events (around midnight) have two dates because casts started pm
% ended am -- use the 1st date.
stations_ = table(unique(stations.Label));
stations_.Properties.VariableNames = {'Label'};
for i = 1:height(stations_)
    d = sort(stations.Date(strcmp(stations.Label, stations_.Label{i})));
    t = sort(stations.Time(strcmp(stations.Label, stations_.Label{i})));
    stations_.Date(i) = d(1);
    stations_.Time(i) = t(1);
end
dat = innerjoin(dat, stations_, 'Keys', 'Label');
% rename columns
dat.Properties.VariableNames({'Longitude_E', 'Latitude_N'}) = ... 
    {'Longitude', 'Latitude'};
% include extra timing columns
dat.t = datenum(dat.Date);
[dat.Year, ~] = datevec(dat.t);
dat.Yearday = yearday(dat.t);
% filter out columns
dat = dat(:,{'Label','t','Date','Time','Year','Yearday','Latitude','Longitude','Depth', ... 
    'NO3_NO2_C','QF_NO3_NO2_C','PO4_C','QF_PO4_C','Si_OH4_C','QF_Si_OH4_C'});
% remove rows flagged as (potentially) poor quality measurements
dat.Properties.VariableNames(contains(dat.Properties.VariableNames, '_C') ...
    & ~contains(dat.Properties.VariableNames, 'QF_')) = {'N','P','Si'};
dat.Properties.VariableNames(contains(dat.Properties.VariableNames, 'QF_')) ...
    = {'Flag_N','Flag_P','Flag_Si'};
dat.N(dat.Flag_N ~= 0) = nan; dat.Flag_N = [];
dat.P(dat.Flag_P ~= 0) = nan; dat.Flag_P = [];
dat.Si(dat.Flag_Si ~= 0) = nan; dat.Flag_Si = [];
% convert from short to long format
dat = stack(dat, {'N','P','Si'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat.Variable = cellstr(dat.Variable);
dat(isnan(dat.Value),:) = [];
dat.Type = repmat({'inorganic'}, [height(dat) 1]);
% merge with PS99 data
dat_nut = [dat_nut; dat];



%% Organic nutrient

%~~~~~~~~~~~~~~~~~~~~~
% PS99 & PS107 cruises
%~~~~~~~~~~~~~~~~~~~~~

OMObsFile = 'copy_data_Engel_etal2019.csv';

dat = readtable(fullfile(obsDir,OMObsFile), 'Format', 'auto');

% Extract the PS99 & PS107 cruise data (omitting everything pre-2016)
dat = dat(contains(dat.EventNumber, {'PS99','PS107'}),:);

% The lat-long columns have inconsistent format -- contain degree symbols 
% and measurements both in decimal-degrees and degree-minutes. Reformat...
ind = ~cellfun('isempty', strfind(dat.Latitude, char(176)));
degr = split(dat.Latitude(ind),char(176)); % convert degree-minutes/seconds into degree-decimal
minute = split(degr(:,2),'.');
degr = cellfun(@(x)str2double(x),degr(:,1));
second = cellfun(@(x)str2double(x),minute(:,2));
minute = cellfun(@(x)str2double(x),minute(:,1));
lat = cellstr(num2str(degr + minute ./ 60 + second ./ 3600));
dat.Latitude(ind) = lat;
dat.Latitude = cellfun(@(x)str2double(x), dat.Latitude);
ind = ~cellfun('isempty', strfind(dat.Longitude, char(176)));
west = contains(dat.Longitude(ind), 'W');
dat.Longitude = strrep(dat.Longitude,'E00','');
dat.Longitude = strrep(dat.Longitude,'E0','');
dat.Longitude = strrep(dat.Longitude,'E','');
dat.Longitude = strrep(dat.Longitude,'W00','');
dat.Longitude = strrep(dat.Longitude,'W0','');
dat.Longitude = strrep(dat.Longitude,'W','');
degr = split(dat.Longitude(ind),char(176));
minute = split(degr(:,2),'.');
degr = cellfun(@(x)str2double(x),degr(:,1));
second = cellfun(@(x)str2double(x),minute(:,2));
minute = cellfun(@(x)str2double(x),minute(:,1));
lon = degr + minute ./ 60 + second ./ 3600;
lon(west) = -lon(west);
dat.Longitude(ind) = cellstr(num2str(lon));
dat.Longitude = cellfun(@(x)str2double(x), dat.Longitude);
% Event format is inconsistent... reformat
dat.Properties.VariableNames{'EventNumber'} = 'Label';
dat.Label = strrep(dat.Label, '/', '_');
dat.Label = strrep(dat.Label, '-', '_');
x = regexp(dat.Label, '_');
for i = 1:length(x)
    dat.Label{i}(x{i}(1)) = '/';
end
dat.Label = strrep(dat.Label, '_', '-');
dat.Label = strrep(dat.Label, ' ', '');
% include leading zeros
x = split(dat.Label, '/');
x2 = x(:,2);
x_ = split(x2, '-');
x_1 = x_(:,1);
x_2 = x_(:,2);
for i = 1:length(x_1)
    n1 = length(x_1{i});
    n2 = length(x_2{i});
    if n1 == 1, x_1{i} = ['00' x_1{i}]; end
    if n1 == 2, x_1{i} = ['0' x_1{i}]; end
    if n2 == 1, x_2{i} = ['0' x_2{i}]; end
end
x_ = join([x_1, x_2], '-');
dat.Label = join([x(:,1), x_], '/');
% include extra timing columns
dat.t = datenum(dat.Date);
dat.Yearday = yearday(dat.t);
% filter out data columns
dat = dat(:,{'Label','t','Date','Year','Yearday','Latitude','Longitude','Depth','chl_a','POC','PON'});
% Convert POM units from (mu g / L) to (m mol / m^3)
dat.PON = dat.PON ./ 14;
dat.POC = dat.POC ./ 12;
% Chl-a units from mu g / L -> m g / m^3, no conversion necessary
% convert from short to long format
dat = stack(dat, {'chl_a','POC','PON'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat.Variable = cellstr(dat.Variable);
dat(isnan(dat.Value),:) = [];
dat.Type = repmat({'organic'}, [height(dat) 1]);

dat_OM = dat;


%~~~~~~~~~~~~~
% PS114 cruise
%~~~~~~~~~~~~~

OMObsFile = 'ARK32-1_PS114_2018_Chla_POC_Nothig_CORRECTED.csv';

dat = readtable(fullfile(obsDir,OMObsFile), 'Format', 'auto');

% reformat lat/long columns -- use +/- insteqad of E/W, and use decimal-degrees
% instead of degree-minutes-seconds
west = contains(dat.Longitude, 'W');
dat.Longitude = strrep(dat.Longitude,'E ','');
dat.Longitude = strrep(dat.Longitude,'W ','');
dat.Longitude = strrep(dat.Longitude,' ','.');
degr = split(dat.Longitude, '.'); % convert degree-minutes/seconds into degree-decimal
minute = str2double(degr(:,2));
second = str2double(degr(:,3));
degr = str2double(degr(:,1));
lon = degr + minute ./ 60 + second ./ 3600;
lon(west) = -lon(west);
dat.Longitude = lon;
dat.Latitude = strrep(dat.Latitude, ' ', '.');
degr = split(dat.Latitude, '.');
minute = str2double(degr(:,2));
second = str2double(degr(:,3));
degr = str2double(degr(:,1));
dat.Latitude = degr + minute ./ 60 + second ./ 3600;
% reformat event labelling...
dat.Properties.VariableNames{'Station'} = 'Label';
x = regexp(dat.Label, '_');
for i = 1:length(x)
    dat.Label{i}(x{i}(1)) = '/';
end
dat.Label = strrep(dat.Label, '_', '-');
dat.Label = strrep(dat.Label, ' ', '');
% include extra timing columns
dat.t = datenum(dat.Date);
[dat.Year,~] = datevec(dat.t);
dat.Yearday = yearday(dat.t);
dat.Properties.VariableNames{'TimeUTC'} = 'Time';
% filter out data columns
dat = dat(:,{'Label','t','Date','Time','Year','Yearday','Latitude','Longitude','Depth','Chla','POC','PON'});
dat.Properties.VariableNames{'Chla'} = 'chl_a';
% Convert POM units from (mu g / L) to (m mol / m^3)
dat.PON = dat.PON ./ 14;
dat.POC = dat.POC ./ 12;
% convert from short to long format
dat = stack(dat, {'chl_a','POC','PON'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat.Variable = cellstr(dat.Variable);
dat(isnan(dat.Value),:) = [];
dat.Type = repmat({'organic'}, [height(dat) 1]);
% merge with PS99 & PS107 data

dat_OM.Time = cell(height(dat_OM),1);
dat_OM = movevars(dat_OM, 'Time', 'After', 'Date');
dat_OM = [dat_OM; dat];


%% Size spectra

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

dat_size.BioVolDensity = 1e-18 * dat_size.cellVolume .* dat_size.cellDensity; % m^3 / m^3

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
    'cellDensity', 'cellDensitySD', 'BioVolDensity'}, strcat('Ndensity_', NsizeEqs), 'Ndensity'];
dat_size = dat_size(:,keepVars);

cols = [[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]]; % plotting colours for different scenarios

v = reshape(varargin, [2 0.5*length(varargin)]);

%~~~
% optional plots
plotNconcSpectra = ~isempty(v) && any(contains(v(1,:),'plotNconcSpectra')) && v{2,strcmp(v(1,:), 'plotNconcSpectra')};
plotCellConcSpectra = ~isempty(v) && any(contains(v(1,:),'plotCellConcSpectra')) && v{2,strcmp(v(1,:), 'plotCellConcSpectra')};
plotBioVolSpectra = ~isempty(v) && any(contains(v(1,:),'plotBioVolSpectra')) && v{2,strcmp(v(1,:), 'plotBioVolSpectra')};

makePlots = table(plotNconcSpectra, plotCellConcSpectra, plotBioVolSpectra);
makePlots.Properties.VariableNames = {'Ndensity', 'cellDensity', 'BioVolDensity'};

if any(makePlots{1,:})
    figure
    nc = sum(makePlots{1,:});
    for i = 1:nc
        pv = makePlots.Properties.VariableNames{:,i};
        subplot(1,nc,i)
        ii = strcmp(dat_size.scenario, scenarios(1));
        loglog(dat_size.ESD(ii), dat_size.(pv)(ii), 'Color', cols(1,:))
        hold on
        for ij = 2:N
            ii = strcmp(dat_size.scenario, scenarios(ij));
            loglog(dat_size.ESD(ii), dat_size.(pv)(ii), 'Color', cols(ij,:))
        end
        gc = gca;
        switch pv
            case 'Ndensity'
                gc.YLim(1) = 1e-6;
                ylab = 'Nitrogen density';
                yunit = 'mmol N m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
            case 'cellDensity'
                gc.YLim(1) = 1e1;
                ylab = 'cell conc. density';
%                 yunit = 'cells m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$';
                yunit = 'cells m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
            case 'BioVolDensity'
                gc.YLim(1) = 1e-3 * 1e-6;
                ylab = 'bio vol density';
                yunit = 'm$^3\,$m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
        end
        gc.YLim(2) = 2 * max(dat_size.(pv));
        xlabel('ESD (\mum)')
        ylabel([ylab ' (' yunit ')'], 'Interpreter', 'latex')
        if i == 1
            for ii = N:-1:1
                j = N-ii+1;
                base = 5; space = 5;
                yt = space^(j-1) * base * gc.YLim(1);
                text(gc.XLim(1) * (1 + 0.3 * 10), yt, scenarios{ii});
                line([gc.XLim(1) * (1 + 0.1 * 10) gc.XLim(1) * (1 + 0.2 * 10)], [yt yt], 'Color', cols(ii,:))
            end
        end
        hold off
    end
    sgtitle('Size spectra data')
end
%~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Merge data sets: group all the 'scalar' data into single table and store
% the size data separately.

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

% Separate cruise, station, and sample from label
x1 = regexp(dat.Label, '/');
x2 = regexp(dat.Label, '-');
for i = 1:length(x1)
    dat.Cruise{i} = dat.Label{i}(1:x1{i}-1);    
    dat.Station{i} = dat.Label{i}(x1{i}+1:x2{i}-1);
    dat.Sample{i} = dat.Label{i}(x2{i}+1:end);
end
dat = movevars(dat, {'Cruise','Station','Sample'}, 'After', 'Label');

% Order rows by time, then station, then sample
x = [dat.t, str2double(string(dat.Station)), str2double(string(dat.Sample))];
[~, datOrd] = sortrows(x);
dat = dat(datOrd,:);

% Index each unique sampling event...
Label = unique(dat.Label, 'stable');
Event = (1:length(Label))';
events = table(Event, Label);
dat = join(dat, events);
dat = movevars(dat, 'Event', 'After', 'Sample');


% Consistency checks...
% Different data sets (inorganic & organic nutrients) sampled from the same
% cruises should have consistent metadata such that measurements (space-time
% coords) from the same sampling events are identical. This is not
% guarenteed because metadata is recorded by different people at slightly
% different times...
cruises = unique(dat.Cruise, 'stable');
for i = 1:length(cruises)    
    di = dat(strcmp(dat.Cruise, cruises{i}),:);
    dataTypes = unique(di.Type);
    if length(dataTypes) > 1
        events = unique(di.Event);
        for j = 1:length(events)
            dj = di(di.Event == events(j),:);
            dataTypes_j = unique(dj.Type);
            if length(dataTypes_j) > 1
                dj.Latitude = repmat(mean(unique(dj.Latitude)), [height(dj) 1]);
                dj.Longitude = repmat(mean(unique(dj.Longitude)), [height(dj) 1]);
                emptyTime = cellfun(@(x) isempty(x), dj.Time);
                if all(~emptyTime)
                    dj.Time = cellstr(repmat(dj.Time{1}, size(dj.Time)));
                end
                if any(emptyTime) && ~all(emptyTime)
                    times = unique(dj.Time(~emptyTime));
                    dj.Time(emptyTime) = repmat(times(1), [sum(emptyTime) 1]);
                end
            end
            di(di.Event == events(j),:) = dj;
        end
    end
    dat(strcmp(dat.Cruise, cruises{i}),:) = di;
end


% Store fitting data as a struct
scalarData = table2struct(dat, 'ToScalar', true);
scalarData.Variable = cellstr(scalarData.Variable);
scalarData.nSamples = size(scalarData.Value,1);
scalarData.nEvents = length(unique(scalarData.Event));

sizeData = table2struct(dat_size, 'ToScalar', true);

Data.scalar = scalarData;
Data.size = sizeData;



end

