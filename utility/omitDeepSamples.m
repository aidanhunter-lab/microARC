function dat = omitDeepSamples(Data, FixedParams, varargin)

% Discard in-situ samples collected from below the maximum modelled depth

dat = Data;
scalarData = Data.scalar;

% remove = scalarData.Depth > max(abs(FixedParams.z)); % omit data sampled from below midpoint of deepest modelled layer
remove = scalarData.Depth > max(abs(FixedParams.zw)); % omit data sampled from below lower bound of deepest modelled layer
if any(remove)
    fields = fieldnames(scalarData);
    for i = 1:length(fields)
        if size(scalarData.(fields{i}), 1) > 1
            scalarData.(fields{i})(remove) = [];
        end
    end
    scalarData.nSamples = length(scalarData.Year);
    scalarData.nEvents = length(unique(scalarData.Event));
end

dat.scalar = scalarData;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Optional arguments to omit samples from below selected depths

extractVarargin(varargin)

if ~exist('chlSampleDepthLimit', 'var')
    chlSampleDepthLimit = inf;
end
if ~exist('DINSampleDepthLimit', 'var')
    DINSampleDepthLimit = inf;
end
if ~exist('PONSampleDepthLimit', 'var')
    PONSampleDepthLimit = inf;
end
if ~exist('POCSampleDepthLimit', 'var')
    POCSampleDepthLimit = inf;
end

ds = dat.scalar;
ds = rmfield(ds, {'nSamples','nEvents','obsInCostFunction'});
ds = struct2table(ds);
omit = (strcmp(ds.Variable, 'chl_a') & ds.Depth >= chlSampleDepthLimit) | ...
    (strcmp(ds.Variable, 'N') & ds.Depth >= DINSampleDepthLimit) | ...
    (strcmp(ds.Variable, 'PON') & ds.Depth >= PONSampleDepthLimit) | ...
    (strcmp(ds.Variable, 'POC') & ds.Depth >= POCSampleDepthLimit);
ds(omit,:) = [];
ds = table2struct(ds, 'ToScalar', true);
ds.nSamples = length(ds.Value);
ds.nEvents = length(unique(ds.Event));
ds.obsInCostFunction = dat.scalar.obsInCostFunction;
dat.scalar = ds;
clear ds

