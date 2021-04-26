function [Data, Forc, FixedParams] = selectYears(Data, Forc, FixedParams, varargin)
% Select which years to model and filter the other years out of the data

extractVarargin(varargin)

% years available for each data type
nut_years = unique(Data.scalar.Year(strcmp(Data.scalar.Type, 'inorganic')));
OM_years = unique(Data.scalar.Year(strcmp(Data.scalar.Type, 'organic')));
size_years = unique(Data.size.Year);
forc_years = Forc.years(:);
all_years = unique([nut_years(:); OM_years(:); size_years(:)]);

% if there is no year where all data types are available then we need to
% aggregate data types across different years -- the last data type 
% specified in prefOrd is most likely to come from different year (forc
% must be in 1st position)
prefOrd = {'forc', 'size', 'OM', 'nut'};
datYear = cell(size(prefOrd)); % datYear stores the year(s) used for each data type

% index which years are present in each data set
mat = false(length(all_years),length(prefOrd)); % [years, data type]
for i = 1:length(prefOrd)
    mat(:,i) = ismember(all_years, eval([prefOrd{i} '_years']));
end

% do any years have all fitting-data types?
anyCompleteYears = any(sum(mat, 2) == length(prefOrd));

if ~exist('singleYear', 'var')
    singleYear = true; % return data from a single year? (default true... false will be trickier and would need a multi-year cost function)
end


mat(~mat(:,1),:) = false; % Forcing data is always required
dataCount = sum(mat(:,2:end), 2);
mostData = dataCount == max(dataCount); % index years with most complete data
% yrs = mostData; % years with most complete data

% Fill datYear wih values dependent upon singleYear and anyCompleteYears
switch singleYear

    case true
    
        switch anyCompleteYears
            
            case true
                % Select the latest year with complete data (rows of 1 in
                % mat) and discard the rest.
                ind = find(mostData, 1, 'last');
                bestYear = all_years(ind);
                for i = 1:length(datYear)
                    datYear{i} = bestYear;
                end
                
            case false                
                % Use the year with the most complete data, then infill
                % remaining missing data from another year.
                tiedYears = sum(mostData) > 1; % do multiple years have equal amounts of data types?

                if tiedYears % then choose between them using prefOrd
                    mat_ = cumsum(mat, 2);
                    [~, I] = max(mat_, [], 2);
                    I(~mostData) = inf;
                    bestYear_ind = I == min(I); % this works because columns of mat are ordered by prefOrd
                    bestYear = all_years(bestYear_ind);
                    missingDataType = prefOrd(~mat(bestYear_ind,:)); % data types not available in bestYear
                    for i = 1:length(datYear)
                        if ~ismember(missingDataType, prefOrd{i})
                            datYear{i} = bestYear;
                        end
                    end
                    for i = 1:length(missingDataType)
                        iy = eval([missingDataType{i}, '_years']); % years available for missing data type
                        dy = abs(iy - bestYear); % proximity to target year
                        dy = dy == min(dy); % choose year closest to target year
                        if sum(dy) > 1
                            % if multiple choices are equally close to target year then choose the later year
                            id = find(dy == 1, 1, 'last');
                            dy = false(size(dy));
                            dy(id) = true;
                        end
                        datYear{strcmp(missingDataType{i}, prefOrd)} = iy(dy);
                    end
                    
                else % no tied years => a bit simpler
                    bestYear_ind = mostData;
                    bestYear = all_years(bestYear_ind);
                    missingDataType = prefOrd(~mat(bestYear_ind,:)); % data types not available in bestYear
                    for i = 1:length(datYear)
                        if ~ismember(missingDataType, prefOrd{i})
                            datYear{i} = bestYear;
                        end
                    end
                    for i = 1:length(missingDataType)
                        iy = eval([missingDataType{i}, '_years']); % years available for missing data type
                        dy = abs(iy - bestYear); % proximity to target year
                        dy = dy == min(dy); % choose year closest to target year
                        if sum(dy) > 1
                            % if multiple choices are equally close to target year then choose the later year
                            id = find(dy == 1, 1, 'last');
                            dy = false(size(dy));
                            dy(id) = true;
                        end
                        datYear{strcmp(missingDataType{i}, prefOrd)} = iy(dy);
                    end
                end
        end
        
    case false
        
        warning('More coding required for multiple target years! More work to do...')
end




% Filter out unused years of forcing data
iy = ismember(forc_years, bestYear);
[Y,~] = datevec(Forc.t);
ind = any(ismember(Y, Forc.years(iy)));
ntraj = length(ind);
ntraj_new = sum(ind);
fields = fieldnames(Forc);
for i = 1:length(fields)
    x = Forc.(fields{i});
    s = size(x);
    if ~any(s == ntraj), continue; else
        nd = length(s);
        ind_ = repmat(reshape(ind, [ones(1,nd-1) ntraj]), [s(1:end-1) 1]);
        Forc.(fields{i}) = reshape(x(ind_), [s(1:end-1), ntraj_new]);
    end
end
Forc.years = Forc.years(iy);
FixedParams.years = FixedParams.years(iy);

% Filter the scalar data
% organic matter
keep = ismember(Data.scalar.Year, datYear{strcmp(prefOrd, 'OM')}) & ...
    strcmp(Data.scalar.Type, 'organic');
% inorganic nutrient
keep = keep | (ismember(Data.scalar.Year, datYear{strcmp(prefOrd, 'nut')}) & ...
    strcmp(Data.scalar.Type, 'inorganic'));

fields = fieldnames(Data.scalar);
for i = 1:length(fields)
    if size(Data.scalar.(fields{i}), 1) ==  Data.scalar.nSamples
        Data.scalar.(fields{i})(~keep) = [];
    end
end
Data.scalar.nSamples = length(Data.scalar.Value);
Data.scalar.nEvents = length(unique(Data.scalar.Event));

% Filter size spectra data
if ~exist('singleSpectra', 'var')
    singleSpectra = true;
end

% data aggregated over events

keep = ismember(Data.size.Year, datYear{strcmp(prefOrd, 'size')});
fields = fieldnames(Data.size);
for i = 1:length(fields)
    if size(Data.size.(fields{i}), 1) == Data.size.nSamples
       Data.size.(fields{i})(~keep) = []; 
    end
end
Data.size.nSamples = size(Data.size.ESD,1);

if singleSpectra
    uscenario = unique(Data.size.scenario);
    if length(uscenario) > 1
        keep = strcmp(Data.size.scenario, uscenario{1}); % this could be better coded, but it works given the data we have!
        for i = 1:length(fields)
            if size(Data.size.(fields{i}), 1) == Data.size.nSamples
                Data.size.(fields{i})(~keep) = [];
            end
        end
    end
end
Data.size.nSamples = size(Data.size.ESD,1);

keep = ismember(Data.size.dataBinned.Year, datYear{strcmp(prefOrd, 'size')});
if singleSpectra
    uscenario = unique(Data.size.dataBinned.scenario(keep));
    if length(uscenario) > 1
        keep = keep & strcmp(Data.size.dataBinned.scenario, uscenario{1});
    end
end
Data.size.dataBinned = structfun(@(x) x(keep), Data.size.dataBinned, ...
    'UniformOutput', false);


% full data separated over events

keep = ismember(Data.sizeFull.Year, datYear{strcmp(prefOrd, 'size')});
fields = fieldnames(Data.sizeFull);
nSamples = size(Data.sizeFull.Year,1);
for i = 1:length(fields)
    if size(Data.sizeFull.(fields{i}), 1) == nSamples
       Data.sizeFull.(fields{i})(~keep,:) = []; 
    end
end
Data.sizeFull.nSamples = size(Data.sizeFull.ESD,1);

if singleSpectra
    uscenario = unique(Data.sizeFull.Cruise);
    if length(uscenario) > 1
        keep = strcmp(Data.sizeFull.scenario, uscenario{1}); % this could be better coded, but it works given the data we have!
        for i = 1:length(fields)
            if size(Data.sizeFull.(fields{i}), 1) == Data.sizeFull.nSamples
                Data.sizeFull.(fields{i})(~keep) = [];
            end
        end
    end
end
Data.sizeFull.nSamples = size(Data.sizeFull.ESD,1);

keep = ismember(Data.sizeFull.dataBinned.Year, datYear{strcmp(prefOrd, 'size')});
if singleSpectra
    uscenario = unique(Data.sizeFull.dataBinned.Cruise(keep));
    if length(uscenario) > 1
        keep = keep & strcmp(Data.sizeFull.dataBinned.Cruise, uscenario{1});
    end
end
Data.sizeFull.dataBinned = structfun(@(x) x(keep), Data.sizeFull.dataBinned, ...
    'UniformOutput', false);


% Relabel the event numbers
eventLabs = table((1:Data.scalar.nEvents)', unique(Data.scalar.Event));
eventLabs.Properties.VariableNames = {'newEvent', 'Event'};
tmp = table(Data.scalar.Event); tmp.Properties.VariableNames = {'Event'};
tmp = join(tmp, eventLabs);
Data.scalar.Event = tmp.newEvent;

tmp = table(Data.sizeFull.Event); tmp.Properties.VariableNames = {'Event'};
tmp = join(tmp, eventLabs);
Data.sizeFull.Event = tmp.newEvent;

tmp = table(Data.sizeFull.dataBinned.Event); tmp.Properties.VariableNames = {'Event'};
tmp = join(tmp, eventLabs);
Data.sizeFull.dataBinned.Event = tmp.newEvent;


