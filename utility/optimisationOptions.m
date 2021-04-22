function [FixedParams, Params, Forc] = optimisationOptions(FixedParams, Params, Forc, varargin)
% Choose tuning parameters and cost function and numerical tuning algorithm
% and any other options related to optimisation can be included here...

extractVarargin(varargin)

% Select which parameters to optimise.
% Choose from the lists: Params.scalars & Params.sizeDependent.
parnames = {'wPOM', 'wp_a', 'wp_b', 'rDON', 'rPON', 'rPOC', ...
    'aP', 'm_b', 'Gmax_a', 'Gmax_b', 'k_G', 'pmax_a', 'pmax_b', ... 
    'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', 'Qmax_delQ_b', ... 
    'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b'};

% Check that all chosen parameters exist -- error if not
if ~all(ismember(parnames, Params.scalars) | ... 
        ismember(parnames, Params.sizeDependent))
    error('Error in "optimisationOptions.m". Invalid choice for "parnames": all tuning parameter names must appear in "Params.scalars" or "Params.sizeDependent".')
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

% Select optimising algorithm (so far only ga is available)
optimiserChoices = {'ga','muga'};
if ~exist('optimiser', 'var')
    optimiser = optimiserChoices{1};
end
FixedParams.optimiser = optimiser;
optimise = str2func(optimiser);
assignin('caller', 'optimise', optimise) % assign optimising algorithm to the workspace

% Select cost function.
% There's a few options for the cost function. Not yet sure which is the
% best... from this list choose either option 3 or 4
costFunctionChoices = { ...
    'LSS', ...
    'RMS', ...
    'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal', ...
    'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet', ...
    'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet' ...
    };

if ~exist('costFunctionType', 'var')
    % costFunctionChoices should be given shorter names...
    costFunctionType = costFunctionChoices{4}; % select cost function
end
FixedParams.costFunction = costFunctionType;
assignin('caller', 'costFunctionLabel', costFunctionType)

% Halt integrations at each trajectory's final sampling events to save time
Forc.integrateFullTrajectory = false;



% Parameters are optimised using a numerical population-based algorithm
if ~exist('popSize', 'var')
    popSize = 100; % number of parameter sets in algorithm population    
end
if ~exist('niter', 'var')
    niter = 10; % algorithm iterations
end
assignin('caller', 'popSize', popSize)
assignin('caller', 'niter', niter)

% Optimising algorithm options
switch optimiser
    case 'ga'
        optimiserOptions = optimoptions('ga', ...
            'PopulationSize', popSize, ...
            'Generations', niter, ...
            'InitialPopulationRange', [lb;ub], ...
            'InitialPopulationMatrix', par0, ...
            'Display', 'iter', ...
            'OutputFcn', @gaStoreHistory, ... % gaStoreHistory.m function dynamically stores output, which is available in workspace even if algorithm is halted prematurely
            'PlotFcn', @gaplotbestf);
end
assignin('caller', 'optimiserOptions', optimiserOptions)















