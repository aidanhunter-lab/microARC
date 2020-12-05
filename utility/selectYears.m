function [Data, Forc, FixedParams] = selectYears(data, forc, fixedpars)
% Select which years to model and filter the other years out of the data
Data = data;
Forc = forc;
FixedParams = fixedpars;

nut_years = unique(Data.scalar.Year(strcmp(Data.scalar.Type, 'inorganic')));
OM_years = unique(Data.scalar.Year(strcmp(Data.scalar.Type, 'organic')));
size_years = unique(Data.size.Year);
forc_years = Forc.years;
all_years = unique([nut_years(:); OM_years(:); size_years(:)]);
% nyears = length(all_years);
% index which years are present in each data set
mat = [ismember(all_years, forc_years), ...
    ismember(all_years, nut_years), ...
    ismember(all_years, OM_years), ...
    ismember(all_years, size_years)];
% Select all years with complete data (rows of 1 in mat) and discard the
% rest. If years with complete data are not present then return the year
% with the most complete data, and infill with other years where needed.
mostData = sum(mat(:,2:end), 2);
mostData = mostData == max(mostData);
yrs = mostData & mat(:,1) == 1;
whichYears = all_years(yrs);

% Filter out unused years of forcing data
iy = ismember(forc_years, whichYears);
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
iy = ismember(OM_years, whichYears);
if ~all(iy == 0)
    keep = ismember(Data.scalar.Year, OM_years(iy)) & strcmp(Data.scalar.Type, 'organic');
else
    keep = strcmp(Data.scalar.Type, 'organic');
end
% inorganic nutrient
iy = ismember(nut_years, whichYears);
if ~all(iy == 0)
    keep = keep | (ismember(Data.scalar.Year, nut_years(iy)) & strcmp(Data.scalar.Type, 'inorganic'));
else
    keep = keep | strcmp(Data.scalar.Type, 'inorganic');
end

fields = fieldnames(Data.scalar);
for i = 1:length(fields)
    if size(Data.scalar.(fields{i}), 1) ==  Data.scalar.nSamples
        Data.scalar.(fields{i})(~keep) = [];
    end
end
Data.scalar.nSamples = length(Data.scalar.Value);
Data.scalar.nEvents = length(unique(Data.scalar.Event));

% Filter size spectra data
iy = ismember(Data.size.Year, whichYears);
fields = fieldnames(Data.size);
for i = 1:length(fields)
    if size(Data.size.(fields{i}), 1) == Data.size.nSamples
       Data.size.(fields{i}) = Data.size.(fields{i})(iy); 
    end
end
Data.size.nSamples = size(Data.size.ESD,1);
iy = ismember(Data.size.dataBinned.Year, whichYears);
Data.size.dataBinned = structfun(@(x) x(iy), Data.size.dataBinned, ...
    'UniformOutput', false);

% Relabel the event numbers
eventLabs = table((1:Data.scalar.nEvents)', unique(Data.scalar.Event));
eventLabs.Properties.VariableNames = {'newEvent', 'Event'};
tmp = table(Data.scalar.Event); tmp.Properties.VariableNames = {'Event'};
tmp = join(tmp, eventLabs);
Data.scalar.Event = tmp.newEvent;

