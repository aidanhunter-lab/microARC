function Data = prepareFittingData(obsDir, varargin)

% Load and clean fitting data, output as struct

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Inorganic nutrient

nutrientObsFile = 'PS99_2_nutrients_rev.csv';

dat_nut = readtable(fullfile(obsDir,nutrientObsFile));
date = split(dat_nut.Date, 'T');
dat_nut.Date = date(:,1);
dat_nut.Time = date(:,2);
dat_nut = movevars(dat_nut, 'Time', 'After', 'Date');
dat_nut.Year = year(dat_nut.Date); % this 'year' function apparently causes problems on some MatLab versions
% dat_nut.Year = cellfun(@(x) str2double(x), ...
%     cellfun(@(x) x(1:4), dat_nut.Date, 'UniformOutput', false));
dat_nut = movevars(dat_nut, 'Year', 'After', 'Time');
dat_nut.Yearday = yearday(datenum(dat_nut.Date));
dat_nut = movevars(dat_nut, 'Yearday', 'After', 'Year');
% filter out data columns
dat_nut = dat_nut(:,{'Year','Yearday','Event','Event_2','Latitude','Longitude','Depth','NO3_NO2_corrected','Flag_NO3_NO2_c'});
% remove rows flagged as (potentially) poor quality measurements
dat_nut(dat_nut.Flag_NO3_NO2_c ~= 0,:) = [];
dat_nut.Flag_NO3_NO2_c = [];
dat_nut.Properties.VariableNames{'NO3_NO2_corrected'} = 'N';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% POM and Chl_a

OMObsFile = 'copy_data_Engel_etal2019.csv';

dat_OM = readtable(fullfile(obsDir,OMObsFile));
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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Size spectra

sizeSpectraObsFile = 'S1-S4_spectra_Nbiomass.csv';

dat_size = readtable(fullfile(obsDir, sizeSpectraObsFile));

NsizeEqs = {'MDL','Moal','Montagnes','Verity'}; % 4 different cell size -> N conversions
ne = length(NsizeEqs);
ESDmin = 2; % min/max sizes to retain from the data
ESDmax = 200;
% xl = [ESDmin ESDmax]; % plot axis limits
yl = [1e-6 1e8];

% Convert cell density to nitrogen density
esdScale = mean(diff(log10(unique(dat_size.ESD)))); % log10 scale of ESD measurement intervals
sf = 1e-6 .* 1/14;  % unit conversion factor: pg/L -> mmol/m^3
dat_size.Ndensity = sf .* esdScale .* dat_size.cellDensity .* ... 
    sum(dat_size{:,strcat('NperCell_', NsizeEqs)},2) ./ ne;
dat_size.NdensitySD = sf .* esdScale ./ ne .* dat_size.cellDensitySD .* ... 
    (sum(dat_size{:,strcat('NperCell_', NsizeEqs)} .^ 2,2)) .^ 0.5;
% dat_size.Ndensity = sf .* dat_size.cellDensity .* ... 
%     sum(dat_size{:,strcat('NperCell_', NsizeEqs)},2) ./ ne;
% dat_size.NdensitySD = sf ./ ne .* dat_size.cellDensitySD .* ... 
%     (sum(dat_size{:,strcat('NperCell_', NsizeEqs)} .^ 2,2)) .^ 0.5;

% Remove unnecessaary variables
% Dat_size = dat_size;
scenarios = unique(dat_size.scenario)'; % data collected from different cruises/years
TrophicLevel = 'autotroph';
dat_size = dat_size(strcmp(dat_size.trophicLevel,TrophicLevel) & ... 
    dat_size.ESD >= ESDmin & dat_size.ESD <= ESDmax, ...
    {'season','ESD','scenario','Ndensity'});
% remove scenario S4 because it's from a winter cruise and the data has a different shape
dat_size(strcmp(dat_size.scenario,scenarios(4)),:) = [];
scenarios(4) = [];
N = length(scenarios);
cols = [[1 0 0]; [0 1 0]; [0 0 1]];

v = reshape(varargin, [2 0.5*length(varargin)]);

if ~isempty(v) && any(contains(v(1,:),'plotRawSizeSpectra')) && v{2,strcmp(v(1,:), 'plotRawSizeSpectra')}
    figure
    ii = strcmp(dat_size.scenario, scenarios(1));
    loglog(dat_size.ESD(ii), dat_size.Ndensity(ii), 'Color', cols(1,:))
    hold on
    for i = 2:N
        ii = strcmp(dat_size.scenario, scenarios(i));
        loglog(dat_size.ESD(ii), dat_size.Ndensity(ii), 'Color', cols(i,:))
    end
    gc = gca; gc.YLim(1) = yl(1); gc = gca;
    xlabel('ESD (\mum)')
    ylabel('Nitrogen density (mmol N m^{-3})')
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


% Partition the size spectra measurements into size-class bins.
% Find the across-scenario (log scale) mean nitrogen spectra (different scenarios may have different numbers of data points...)
allESD = unique(dat_size.ESD);
allESDn = length(allESD);
allESD = table((1:allESDn)', allESD); allESD.Properties.VariableNames = {'index','ESD'};
mat = nan(height(allESD),N);
for i = 1:N
    ind = strcmp(dat_size.scenario, scenarios(i));
    dt = dat_size(ind, {'ESD','Ndensity'});    
    dtt = innerjoin(allESD, dt);
    mat(:,i) = dtt.Ndensity;
end
dt = table(repmat({'null'}, [allESDn 1]), allESD.ESD, repmat({'mean'}, [allESDn 1]), ...
    10 .^ nanmean(log10(mat), 2));
dt.Properties.VariableNames = dat_size.Properties.VariableNames;
dat_size = [dat_size; dt];

ind_mean = strcmp(dat_size.scenario,'mean');
ibin = partitionSizeSpectra(dat_size.Ndensity(ind_mean));

plotSizeClassBins = 0;
if plotSizeClassBins
    figure
    loglog(dat_size.ESD(ind_mean), dat_size.Ndensity(ind_mean), 'Color', [0 0 0])
    hold on
    xlabel('ESD (\mum)')
    ylabel('Nitrgen (mmol N m^{-3})')
%     ylabel('Nitrgen (\mug L^{-1})')
    title({'mean size spectra:', 'intervals bounded by turning points'})
    gc = gca; yt = gc.YLim;
    for i = 1:length(ibin)
        xt = allESD.ESD(ibin(i));
        line([xt xt], yt, 'Color', [0 0 0], 'LineStyle', ':')
    end    
    hold off
end



% Reduce number of size-class bins by incrementally removing the narrowest bins
nbins = length(ibin) - 1;
bins = cell(nbins - 1, 1);
bins{1} = ibin;
nbd = diff(ibin); % number of data points within each bin
[~, o] = sort(nbd);
while nbins > 2
    rb = o(1); % narrowest bin -- remove 1 of its boundaries
    if rb == 1, removeBin = 2; end
    if rb == nbins, removeBin = nbins - 1; end
    if rb > 1 && rb < nbins
        trb = [nbd(rb - 1) nbd(rb + 1)];
        removeBin = rb + (find(trb == min(trb)) - 1);
    end
    ibin(removeBin) = [];
    nbins = length(ibin) - 1;
    bins{length(bins)-nbins+2} = ibin;
    nbd = diff(ibin);
    [~, o] = sort(nbd);
end


plotSizeClassBins = 0;
if plotSizeClassBins
    figure
    bp = length(bins);
    nr = ceil(sqrt(bp));
    nc = ceil(bp / nr);
    ind_mean = strcmp(dat_size.scenario, 'mean');
    for i = 1:bp
        subplot(nr, nc, i)
        loglog(allESD.ESD, dat_size.Ndensity(ind_mean), 'Color', [0 0 0])
        hold on
        vl = bins{i};
        nvl = length(vl);
        xlabel('ESD (\mum)')
        ylabel('Nitrogen (mmol N m^{-3})')
        title([num2str(nvl-1) ' intervals'])
        gc = gca; yl = gc.YLim;
        for j = 1:length(vl)
            line([allESD.ESD(vl(j)) allESD.ESD(vl(j))], yl, 'LineStyle', ':', 'Color', [0 0 0])
        end
        hold off
    end
end

% Choose a number of intervals, by eye for now, and plot within-interval
% distributions. (Then repeat for all intervals to choose between them...)
bns = bins{6};
nb = length(bns) - 1; % number of size-class bins

if ~isempty(v) && any(contains(v(1,:),'plotSizeClassBins')) && v{2,strcmp(v(1,:), 'plotSizeClassBins')}
    figure
    ind_mean = strcmp(dat_size.scenario, 'mean');
    loglog(allESD.ESD, dat_size.Ndensity(ind_mean), 'Color', [0 0 0])
    hold on
    gc = gca; yl = gc.YLim;
    for i = 1:nb+1
        line([allESD.ESD(bns(i)) allESD.ESD(bns(i))], yl, 'Color', [0 0 0], 'LineStyle', ':')
    end
    xlabel('ESD (\mum)')
    ylabel('Nitrogen (mmol N m^{-3})')
    title({'mean size spectra', 'intervals bounded by selected turning points'})
    hold off
end

% Remove the across-scenario mean data
dat_size(strcmp(dat_size.scenario,'mean'),:) = [];
% Group the data by size-class bin
dat_size.sizeClass = nan(height(dat_size), 1);
for i = 1:nb
    if i ~= nb
        ind = allESD.ESD(bns(i)) <= dat_size.ESD & allESD.ESD(bns(i+1)) > dat_size.ESD;
    else
        ind = allESD.ESD(bns(i)) <= dat_size.ESD & allESD.ESD(bns(i+1)) >= dat_size.ESD;
    end
    dat_size.sizeClass(ind) = i;
end

Sizes = allESD.ESD(bns(1:nb)) + 0.5 * diff(allESD.ESD(bns)); % midpoints within each interval
Sizes = round(Sizes * 2) / 2; % round to 0.5
K = length(Sizes);
dat_size.size = Sizes(dat_size.sizeClass);
dat_size.ndat = nan(height(dat_size), 1); % number of data points per size-class bin
for i = 1:N
    ind0 = strcmp(dat_size.scenario, scenarios{i});
    for j = 1:K
        ind = ind0 & dat_size.sizeClass == j;
        ndat = sum(ind);
        dat_size.ndat(ind) = repmat(ndat, [ndat 1]);
    end
end

if ~isempty(v) && any(contains(v(1,:),'plotRawSpectraAndBins')) && v{2,strcmp(v(1,:), 'plotRawSpectraAndBins')}
    figure
    sc = strcmp(dat_size.scenario, scenarios(1));
    loglog(dat_size.ESD(sc), dat_size.Ndensity(sc), 'Color', cols(1,:))
    hold on
    for i = 2:N
        sc = strcmp(dat_size.scenario, scenarios(i));
        loglog(dat_size.ESD(sc), dat_size.Ndensity(sc), 'Color', cols(i,:))        
    end
    gc = gca; 
    if gc.YLim(1) < 1e-6, gc.YLim(1) = 1e-6; end
    yl = gc.YLim;    
    for i = 1:nb+1
        line([allESD.ESD(bns(i)) allESD.ESD(bns(i))], yl, 'Color', [0 0 0], 'LineStyle', ':')
    end
    for i = N:-1:1
        j = N-i+1;        
        base = 5; space = 5;        
        yt = space^(j-1) .* base .* gc.YLim(1);
%         yt = 1*10^j * gc.YLim(1);
        text(gc.XLim(1) * (1 + 0.3 * 10), yt, scenarios{i});
        line([gc.XLim(1) * (1 + 0.1 * 10) gc.XLim(1) * (1 + 0.2 * 10)], [yt yt], 'Color', cols(i,:))
    end
    xlabel('ESD (\mum)')
    ylabel('Nitrogen (mmol N m^{-3})')
    title({'size spectra data', 'intervals bounded by selected turning points'})
    hold off
end

% Model output will comprise nitrogen density estimates for each size-class
% bin, which should be equivalent to summing the data within each bin.
% Thus we can simply multiply the data within each bin by the number of data
% points in that bin, to get data distributions for each bin.
% These can then be scaled for use in a cost function.

dat_size.NdensityTot = nan(height(dat_size), 1);
for i = 1:N
    ind0 = strcmp(dat_size.scenario, scenarios(i));
    for j = 1:nb
        ind = ind0 & dat_size.sizeClass == j;
        n = unique(dat_size.ndat(ind));
        dat_size.NdensityTot(ind) = n .* dat_size.Ndensity(ind);
    end
end


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
scenarioTimes.S1.date = datetime(2016,7,16);
scenarioTimes.S2.date = datetime(2017,8,19);
scenarioTimes.S3.date = datetime(2018,8,3);
scenarioTimes.S1.depth = [10 30];
scenarioTimes.S2.depth = [10 40];
scenarioTimes.S3.depth = [5 43];
dat_size.Year = nan(height(dat_size), 1);
dat_size.Yearday = nan(height(dat_size), 1);
dat_size.DepthMin = nan(height(dat_size), 1);
dat_size.DepthMax = nan(height(dat_size), 1);
fields = fieldnames(scenarioTimes);
for i = 1:length(fields)
    ind = strcmp(dat_size.scenario,fields{i});
    indn = sum(ind);
    [yr, ~] = datevec(scenarioTimes.(fields{i}).date);
    yrDay = yearday(datenum(scenarioTimes.(fields{i}).date));
    dat_size.Year(ind) = repmat(yr, [indn 1]);
    dat_size.Yearday(ind) = repmat(yrDay, [indn 1]);    
    dat_size.DepthMin(ind) = scenarioTimes.(fields{i}).depth(1);
    dat_size.DepthMax(ind) = scenarioTimes.(fields{i}).depth(2);    
end

% dat_size = removevars(dat_size, {'season', 'scenario', 'ESD', 'sizeClass'});
% dat_size = removevars(dat_size, {'season', 'scenario', 'ESD', 'sizeClass', 'Ndensity', 'logNdensityTot', 'mu', 'sig', 'resid'});
% dat_size = removevars(dat_size, {'season', 'scenario', 'ESD', 'sizeClass', 'Ndensity'});
dat_size = movevars(dat_size, {'Year', 'Yearday'}, 'Before', 'season');
dat_size = movevars(dat_size, {'scenario'}, 'Before', 'Year');
dat_size = movevars(dat_size, {'DepthMin', 'DepthMax'}, 'After', 'season');
dat_size = movevars(dat_size, {'Ndensity'}, 'Before', 'NdensityTot');
% dat_size = movevars(dat_size, {'Year', 'Yearday', 'DepthMin', 'DepthMax'}, 'Before', 'size');
dat_size.Properties.VariableNames{'NdensityTot'} = 'N_at_size'; % rename the size data measurement variable

% Convert tables to long-form in the measurement variables, then merge
dat_OM = stack(dat_OM, {'chl_a','POC','PON'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat_nut = stack(dat_nut, 'N', ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat_size = stack(dat_size, 'N_at_size', ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat = [dat_nut; dat_OM];

% MOVE THIS OUTSIDE THE FUNCTION
% % Remove samples from below the maximum modelled depth
% dat = dat(dat.Depth < -min(FixedParams.z),:);
% dat_size = dat_size(dat_size.DepthMin < -min(FixedParams.z),:);

% Remove NaNs & rows with non-positive measures
dat(dat.Value <= 0 | isnan(dat.Value),:) = [];
dat_size(dat_size.Value <= 0 | isnan(dat_size.Value),:) = [];

% Index each unique sampling event
EventLabel = unique(dat.EventLabel, 'stable');
Event = (1:length(EventLabel))';
events = table(Event, EventLabel);
dat = join(dat, events);

% Store fitting data as a struct
scalarData = table2struct(dat, 'ToScalar', true);
scalarData.Variable = cellstr(scalarData.Variable);
scalarData.nSamples = size(scalarData.Value,1);
scalarData.nEvents = length(unique(scalarData.Event));
sizeData = table2struct(dat_size, 'ToScalar', true);
sizeData.Variable = cellstr(sizeData.Variable);
sizeData.nSamples = size(sizeData.Value,1);

Data.scalar = scalarData;
Data.size = sizeData;



end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%~~~~~~~~~~
% Functions
%~~~~~~~~~~

% function y = nanunique(x)
%   y = unique(x);
%   if any(isnan(y))
%     y(isnan(y)) = []; % remove all nans
%   end
% end

function y = partitionSizeSpectra(x)
    % uses stationary points to partition size spectra vectors into size-class
    % intervals
    n = length(x);
    xdim = size(x);
    if length(xdim) ~= 2 || ~any(xdim == 1)
        warning('partitionSizeSpectra only accepts vector arguments')
    else
        y = diff(x(:));
        y = [y(1); y];
        y = y > 0;
        y = (y(1:n-1) & ~y(2:n)) | (~y(1:n-1) & y(2:n));
        y = [0; y];
        y = [1; find(y); n]; % endpoints and turning points
        if xdim(1) == 1
            y = y';
        end
    end
end

% function y = movAv(x,n)
%     % moving average to smooth data to reduce number of size-class 
%     % intervals found with partitionSizeSpectra.
%     % x = data to smooth
%     % n = number of points either side of focal data point
%     N = length(x);
%     y = nan(N,1);
%     y(1) = x(1); y(N) = x(N);
%     span = 2 * n + 1;
%     Ind = (n+1):(N-n);
%     mat = nan(length(Ind),span);
%     v = -n:n;
%     for i = 1:span
%         mat(:,i) = x(Ind + v(i));
%     end
%     av = sum(mat, 2) ./ span;
%     y(Ind) = av;
%     % edge points
%     for i = n:-1:2
%         nn = i - 1;
%         v = -nn:nn;
%         span = 2 * nn + 1;
%         av = sum(x(i + v)) ./ span;
%         y(i) = av;
%     end
%     j = (N-n+1):(N-1);
%     for i = 1:n-1
%         nn = n - i;
%         v = -nn:nn;
%         span = 2 * nn + 1;
%         av = sum(x(j(i) + v)) / span;
%         y(j(i)) = av;
%     end
% end
% 
% function y = rnormTruncated(n,mu,sig,a,b)
%     pa = normcdf(a,mu,sig);
%     y = norminv(pa + rand(n,1) .* (normcdf(b,mu,sig) - pa)) .* sig + mu;
% end


