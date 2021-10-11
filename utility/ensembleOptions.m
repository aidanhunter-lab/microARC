function [FixedParams, Params, Forc, Data] = ... 
    ensembleOptions(FixedParams, Params, Forc, Data, varargin)

% This function is modified from optimisationOptions.m

fitToFullSizeSpectra = [];

extractVarargin(varargin)

if isempty(fitToFullSizeSpectra) || ~islogical(fitToFullSizeSpectra)
    % if unspecified then by default fit model using binned size spectra data
    fitToFullSizeSpectra = false;
end

FixedParams.fitToFullSizeSpectra = fitToFullSizeSpectra;

%% Generate values for parameters that vary in ensemble run

if ~exist('parSuite', 'var')
    error('A cell array, "parSuite", containing parameter names must be included within varargin as a name-value pair. (Technically, this should not feature as an "optional" argument!)')
end
if ~exist('parRes', 'var')
    parRes = 8; % By default generate 8 values for each parameter
end

npars = length(parSuite); % Number of varied parameters
% Number of unique parameter combinations (= number of function evalutions)
nCombinations = parRes .^ npars;

parSuite

% Values generated using latin hypercube
parVals = lhsdesign(parRes, npars);
% Rearrange into npars-dimensional form, with parRes values per dimension






% Choose from the lists: Params.scalars & Params.sizeDependent.
parnames = {'wPOM1', 'wp_a', 'wp_b', 'rDON', 'rPON', ...
    'aP', 'm_a', 'm2', 'Gmax_a', 'Gmax_b', 'k_G', 'pmax_a', 'pmax_b', ... 
    'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', 'Qmax_delQ_b', ... 
    'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b'};

% Check that all chosen parameters exist -- error if not
if ~all(ismember(parnames, Params.scalarParams) | ... 
        ismember(parnames, Params.vectorParams))
    error('Error in "optimisationOptions.m". Invalid choice for "parnames": all tuning parameter names must appear in "Params.scalarParams" or "Params.vectorParams".')
end

% Extract tuning parameter initial values and bounds
npars = length(parnames);
par0 = cell2mat(cellfun(@(x) Params.(x), parnames, 'UniformOutput', false));
bounds = cellfun(@(x) Params.bounds.(x), parnames, 'UniformOutput', false);
lb = cellfun(@(x) x(1), bounds); % lower bounds
ub = cellfun(@(x) x(2), bounds); % upper bounds
% Store tuning parameter names and their bounding values
FixedParams.tunePars = parnames;
FixedParams.tunePars_lb = lb;
FixedParams.tunePars_ub = ub;
% Assign to workspace
assignin('caller', 'npars', npars)
assignin('caller', 'boundsLower', lb)
assignin('caller', 'boundsUpper', ub)

%% Select cost function.
% There's a few options for the cost function. Not yet sure which is the
% best... Hellinger2_groupWaterOrigin is most defensible as it makes fewest
% assumptions
[~, ~, costFunctionChoices] = costFunction();

if ~exist('costFunctionType', 'var')
    costFunctionType = 'meanCDFdist_Hellinger';
    if ~ismember(costFunctionType, costFunctionChoices)
        costFunctionType = costFunctionChoices{1};
    end
end
FixedParams.costFunction = costFunctionType;
assignin('caller', 'costFunctionLabel', costFunctionType)


%% Filter forcing- and fitting-data

if ~exist('fitTrajectories', 'var')
    fitTrajectories = [];
end

[Forc, Data] = filterInputByOrigin(Forc, Data, 'fitTrajectories', fitTrajectories);

