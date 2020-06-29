function Data = prepareFittingData(obsDir, FixedParams)

% Load and clean fitting data, output as struct

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Inorganic nutrient

nutrientObsFile = 'PS99_2_nutrients_rev.csv';

dat_nut = readtable(fullfile(obsDir,nutrientObsFile));
date = split(dat_nut.Date, 'T');
dat_nut.Date = date(:,1);
dat_nut.Time = date(:,2);
dat_nut = movevars(dat_nut, 'Time', 'After', 'Date');
dat_nut.Year = year(dat_nut.Date);
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
% Convert POM units from (mu g / L) to (mu mol / L)
dat_OM.PON = dat_OM.PON ./ 14;
dat_OM.POC = dat_OM.POC ./ 12;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Merge data sets

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

% Convert tables to long-form then merge
dat_OM = stack(dat_OM, {'chl_a','POC','PON'}, ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');
dat_nut = stack(dat_nut, 'N', ... 
    'IndexVariableName','Variable','NewDataVariableName','Value');

dat = [dat_nut; dat_OM];

% Remove samples from below the maximum modelled depth
dat = dat(dat.Depth < -min(FixedParams.z),:);

% Remove NaNs & rows with non-positive measures -- maybe the 'flag' columns
% should be used for filtering...
dat(dat.Value <= 0 | isnan(dat.Value),:) = [];

% Index each unique sampling event
EventLabel = unique(dat.EventLabel, 'stable');
Event = (1:length(EventLabel))';
events = table(Event, EventLabel);
dat = join(dat, events);
dat = movevars(dat, 'Event', 'Before', 'EventLabel');

% Store fitting data as a struct
Data = table2struct(dat, 'ToScalar',true);
Data.Variable = cellstr(Data.Variable);
Data.nSamples = size(Data.Value,1);
Data.nEvents = length(unique(Data.Event));

