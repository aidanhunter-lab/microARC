function Data = prepareFittingData(obsDir, varargin)
% Load and clean fitting data, output as struct

extractVarargin(varargin)

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
% There is no date provided in dat, so need to combine table of sampling 
% station data with table of measurements, but the formating is inconsistent...
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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% size spectra aggregated over sampling events
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% load data
sizeSpectraObsFile = 'S1-S4_spectra_Nbiomass.csv';
dat_size = readtable(fullfile(obsDir, sizeSpectraObsFile), 'Format', 'auto');

% A selection of cell size -> N conversions are contained in these data.
% Although the different conversions produce nitrogen size-spectra of
% similar shape, there is some substantial variation in magnitudes
% (especially at small and large cell sizes).
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
dat_size.cellDensity = 1e3 * dat_size.cellDensity; % cells/L/log10(ESD) -> cells/m^3/log10(ESD)
dat_size.cellDensitySD = 1e3 * dat_size.cellDensitySD; % cells/L/log10(ESD) -> cells/m^3/log10(ESD)
dat_size.cellDensitySE = 1e3 * dat_size.cellDensitySE; % cells/L/log10(ESD) -> cells/m^3/log10(ESD)
dat_size.cellDensityCImin = 1e3 * dat_size.cellDensityCImin; % cells/L/log10(ESD) -> cells/m^3/log10(ESD)
dat_size.cellDensityCImax = 1e3 * dat_size.cellDensityCImax; % cells/L/log10(ESD) -> cells/m^3/log10(ESD)
for i = 1:ne
    % nitrogen per cell: pg N / cell -> mmol N / cell
    dat_size.(['NperCell_' NsizeEqs{i}]) = ... 
        (1/14) * 1e-9 * dat_size.(['NperCell_' NsizeEqs{i}]);
end

% Derive bio-volume variable
scaleFactor = 1; % may alter biovolume density units using scaleFactor... scaleFactor = 1 => units = (mu m)^3 / m^3 /log10(ESD)
dat_size.BioVolDensity = scaleFactor .* dat_size.cellVolume .* dat_size.cellDensity;
dat_size.BioVolDensitySD = abs(scaleFactor .* dat_size.cellVolume) .* dat_size.cellDensitySD;

for i = 1:ne
    % nitrogen density: mug N / L -> mmol N / m^3
    dat_size.(['Ndensity_' NsizeEqs{i}]) = ... 
        (1/14) * dat_size.(['Ndensity_' NsizeEqs{i}]);
end

% Take average values over the various size -> N conversions
dat_size.NperCell = mean(dat_size{:,strcat('NperCell_', NsizeEqs)}, 2);
dat_size.Ndensity = mean(dat_size{:,strcat('Ndensity_', NsizeEqs)}, 2);

% Remove unnecessary variables
% scenarios = unique(dat_size.scenario)'; % data collected from different cruises/years
keepVars = [{'scenario', 'season', 'regime', 'trophicLevel', 'ESD', 'cellVolume', ...
    'cellDensity', 'cellDensitySD', 'BioVolDensity', 'BioVolDensitySD'}, strcat('Ndensity_', NsizeEqs), 'Ndensity'];
dat_size = dat_size(:,keepVars);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% size spectra data measured separately for each event
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% load data
sizeSpectraObsFile = 'sample_spectra_ah.csv';
dat_size_all = readtable(fullfile(obsDir, sizeSpectraObsFile), 'Format', 'auto');

% Convert units of original data into units used in model (keep all
% elemental concentrations in mmol)
dat_size_all.cellDensity = 1e3 * dat_size_all.cellDensity; % cells/L/log10(ESD) -> cells/m^3/log10(ESD)
dat_size_all.BioVolDensity = scaleFactor * dat_size_all.cellVol .* dat_size_all.cellDensity; % (mu m)^3 / m^3 / log10(ESD) if scaleFactor = 1

% Reformat labels
dat_size_all.Properties.VariableNames({'PangaeaEventLabel','depth','lat','long','date'}) = ...
    {'Label','Depth','Latitude','Longitude','Date'};
dat_size_all.Label = strrep(dat_size_all.Label, '_', '/');

% include leading zeros
x = split(dat_size_all.Label, '/');
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
dat_size_all.Label = join([x(:,1), x_], '/');

% include extra timing columns
x = split(dat_size_all.Date, ' ');
dat_size_all.Date = x(:,1);
dat_size_all.t = datenum(dat_size_all.Date);
[dat_size_all.Year,~] = datevec(dat_size_all.t);
dat_size_all.Yearday = yearday(dat_size_all.t);


cruises = unique(dat_size_all.cruise, 'stable'); % data collected from different cruises/years

nutrition = unique(dat_size_all.nutrition, 'stable');


%~~~
% optional plots
if ~exist('plotCellConcSpectra', 'var')
    plotCellConcSpectra = false;
end
if ~exist('plotBioVolSpectra', 'var')
    plotBioVolSpectra = false;
end
if ~exist('plotNconcSpectra', 'var')
    plotNconcSpectra = false;
end

makePlots = table();
switch plotNconcSpectra, case true ,makePlots.Ndensity = true; end
switch plotCellConcSpectra, case true, makePlots.cellDensity = true; end
switch plotBioVolSpectra, case true, makePlots.BioVolDensity = true; end

cluster = unique(dat_size_all.cluster);
cols = [[0 0 1]; [1 0 0]]; % Colour plot lines by water temperature
alpha = 0.35;

if ~isempty(makePlots)
    np = sum(makePlots{1,:}); % number of plots
    nc = length(cruises); % number of columns
    nr = length(nutrition); % number of rows
    for ip = 1:np
        figure
        pv = makePlots.Properties.VariableNames{:,ip};
        xplt = sort(unique(dat_size_all.ESD));
        for ir = 1:nr
            for ic = 1:nc
                subplot(nr, nc, (ir-1)*nc + ic)
                j = strcmp(dat_size_all.cruise, cruises{ic}) & ... 
                    strcmp(dat_size_all.nutrition, nutrition{ir});
                d = dat_size_all(j,:);
                ne = unique(d.Sample, 'stable');
                Ne = length(ne);
                yplt = nan(length(xplt),Ne);
                waterTemp = cell(Ne,1);
                Colour = nan(Ne, 3);
                for ie = 1:Ne
                    jj = d.Sample == ne(ie);
                    waterTemp(ie) = unique(d.cluster(jj));                    
                    Colour(ie,:) = cols(strcmp(waterTemp{ie}, cluster),:);
                    yplt(:,ie) = d.(pv)(jj);                    
                    loglog(xplt, yplt(:,ie), 'Color', [Colour(ie,:) alpha])
                    if ie == 1, hold on; end
                    if ie == ne, hold off; end
                end
                gc = gca;
                switch pv
                    case 'Ndensity'
                        gc.YLim(1) = 1e-6;
                        ylab = 'Nitrogen density';
                        yunit = 'mmol N m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
                        Title = 'Nitrogen conc. density spectra: all sampling events';
                    case 'cellDensity'
                        gc.YLim(1) = 1e1;
                        ylab = 'cell conc. density';
                        %                 yunit = 'cells m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$';
                        yunit = 'cells m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
                        Title = 'Cell conc. density spectra: all sampling events';
                    case 'BioVolDensity'
                        gc.YLim(1) = min(d.cellVol);
                        ylab = 'bio vol density';
                        yunit = '$\mu$m$^3\,$m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
                        Title = 'Bio-volume density spectra: all sampling events';
                end
                gc.YLim(2) = 2 * max(dat_size_all.(pv));
                xlabel('ESD (\mum)')
                if ic == 1
                    ylabel({nutrition{ir}, [ylab ' (' yunit ')']}  , 'Interpreter', 'latex')
                end
                
                if ir == 1
                    title([cruises{ic}, ' (', num2str(unique(dat_size_all.Year(strcmp(dat_size_all.cruise, cruises{ic})))), ')'])
                end
            end
        end
        sgtitle(Title)
    end
    
    % plot the averages over warm/cold regimes
    
    Colour = nan(length(cluster), 3);
    Colour(strcmp(cluster, 'warm'),:) = [1 0 0];
    Colour(strcmp(cluster, 'cold'),:) = [0 0 1];
    
    for ip = 1:np
        figure
        pv = makePlots.Properties.VariableNames{:,ip};
        xplt = sort(unique(dat_size_all.ESD));
        for ir = 1:nr
            for ic = 1:nc
                subplot(nr, nc, (ir-1)*nc + ic)
                j = strcmp(dat_size_all.cruise, cruises{ic}) & ...
                    strcmp(dat_size_all.nutrition, nutrition{ir});
                d = dat_size_all(j,:);
                ne = unique(d.Sample, 'stable');
                Ne = length(ne);
                yplt = nan(length(xplt),Ne);
                waterTemp = cell(Ne,1);
                
                for ie = 1:Ne
                    jj = d.Sample == ne(ie);
                    waterTemp(ie) = unique(d.cluster(jj));
                    yplt(:,ie) = d.(pv)(jj);
                end
                yplt_mean = nan(size(yplt,1), length(cluster));
                for im = 1:length(cluster)
                    ind = strcmp(waterTemp, cluster{im});
                    yplt_mean(:,im) = mean(yplt(:,ind),2);
                end

                h = loglog(xplt, yplt_mean);                
                set(h, {'color'}, num2cell(Colour, 2))

                gc = gca;
                switch pv
                    case 'Ndensity'
                        gc.YLim(1) = 1e-6;
                        ylab = 'Nitrogen density';
                        yunit = 'mmol N m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
                        Title = 'Nitrogen conc. density spectra: warm/cold regime averages';
                    case 'cellDensity'
                        gc.YLim(1) = 1e1;
                        ylab = 'cell conc. density';
                        %                 yunit = 'cells m$^{-3}\,\log_{10}(\frac{ESD}{1\mu m})^{-1}$';
                        yunit = 'cells m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
                        Title = 'Cell conc. density spectra:  warm/cold regime averages';
                    case 'BioVolDensity'
                        gc.YLim(1) = min(d.cellVol);
%                         gc.YLim(1) = 1e-3 * 1e-6;
                        ylab = 'bio vol density';
                        yunit = '$\mu$m$^3\,$m$^{-3}\,\log_{10}($ESD$/1\mu$m$)^{-1}$';
                        Title = 'Bio-volume density spectra:  warm/cold regime averages';
                end
                gc.YLim(2) = 2 * max(dat_size_all.(pv));
                xlabel('ESD (\mum)')
                if ic == 1
                    ylabel({nutrition{ir}, [ylab ' (' yunit ')']}  , 'Interpreter', 'latex')
                end
                
                if ir == 1
                    title([cruises{ic}, ' (', num2str(unique(dat_size_all.Year(strcmp(dat_size_all.cruise, cruises{ic})))), ')'])
                end
                
                if ir == 1 && ic ==1
                    lgd = legend(cluster, 'Location', 'south');
                    title(lgd, 'Regime')
                end

                
            end
        end
        sgtitle(Title)
    end
    
end

%~~~



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Merge data sets: group all the 'scalar' data into single table and store
% the size data separately.

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



% Organise the aggregated size data

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

sizeData = table2struct(dat_size, 'ToScalar', true);



% Organise the full size data table
dat_size_all.Properties.VariableNames{'Sample'} = 'SampleOrig';
% Separate cruise, station, and sample from label
x = split(dat_size_all.Label, '/');
x_ = split(x(:,2), '-');
dat_size_all.Cruise = x(:,1);
dat_size_all.Station = x_(:,1);
dat_size_all.Sample = x_(:,2);
dat_size_all = movevars(dat_size_all, {'Cruise','Station','Sample'}, 'After', 'Label');
dat_size_all = movevars(dat_size_all, {'HGStation'}, 'Before', 'Station');
dat_size_all.cruise = [];
dat_size_all.SampleOrig = [];

% Reformat dates
dat_size_all.Date = datestr(dat_size_all.Date, 'yyyy-mm-dd');
dat_size_all = movevars(dat_size_all, {'t','Date','Year','Yearday'}, 'After', 'Sample');
dat_size_all.SamplingDate = [];

dat_size_all = movevars(dat_size_all, {'Latitude','Longitude'}, 'After', 'Yearday');
dat_size_all.Properties.VariableNames({'cluster','nutrition','cellVol'}) = ...
    {'regime','trophicLevel','cellVolume'};

% Index each unique sampling event... indices should be consistent with the
% scalar data
Label_scalar = unique(dat.Label, 'stable');

% Discard size data with Label not present in the scalar data.
% This is all MSM77 records and 4 stations from PS107.
ind = ismember(dat_size_all.Label, Label_scalar);
dat_size_all = dat_size_all(ind,:);

dat_size_all.order = (1:height(dat_size_all))';

tt = unique(table(dat.Label, dat.Event), 'stable');
tt.Properties.VariableNames = {'Label','Event'};

dat_size_all = innerjoin(dat_size_all, tt);
[~,I] = sort(dat_size_all.order);
dat_size_all = dat_size_all(I,:);

dat_size_all = movevars(dat_size_all, 'Event', 'After', 'Sample');
dat_size_all.order = [];


sizeDataAll = table2struct(dat_size_all, 'ToScalar', true);

Data.scalar = scalarData;
Data.size = sizeData;
Data.sizeFull = sizeDataAll;



end

