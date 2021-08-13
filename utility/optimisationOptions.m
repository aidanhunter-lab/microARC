function [FixedParams, Params, Forc, Data] = ... 
    optimisationOptions(FixedParams, Params, Forc, Data, varargin)
% Choose tuning parameters and cost function and numerical tuning algorithm
% and any other options related to optimisation can be included here...

extractVarargin(varargin)

if ~exist('fitToFullSizeSpectra', 'var') || isempty(fitToFullSizeSpectra)
    % If unspecified then by default fit model using binned size spectra
    % data. Fitting to the full (unbinned) size data creates a very rough
    % cost function surface that's problematic to fit.
    fitToFullSizeSpectra = false;
end
FixedParams.fitToFullSizeSpectra = fitToFullSizeSpectra;

if ~exist('rescaleForOptim', 'var')
    % True by default because I think it's useful to estimate some of the
    % size-dependency (exponent) parameters on a transformed scale.
    % If further transforms were introduced (maybe useful for regularising 
    % parameter search space) then the 'rescaleForOptim' argument will 
    % probably need to be a character string rather than logical.
    rescaleForOptim = true;
end



%% Select parameters to optimise

% Choose from the lists: Params.scalars & Params.sizeDependent.
% parnames = {'wPOM1', 'wp_a', 'wp_b', 'rDON', 'rPON', ...
%     'aP', 'm_a', 'm2', 'Gmax_a', 'Gmax_b', 'k_G', 'pmax_a', 'pmax_b', ... 
%     'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', 'Qmax_delQ_b', ... 
%     'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b'};
parnames = {'wPOM1', 'rDON', 'rPON', 'aP', 'm_a', 'Gmax_a', 'Gmax_b', ... 
    'k_G', 'pmax_a', 'pmax_b', 'Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a', ...
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
% % Assign to workspace
% assignin('caller', 'npars', npars)
% assignin('caller', 'boundsLower', lb)
% assignin('caller', 'boundsUpper', ub)

% Parameter transforms for optimisation search space
switch rescaleForOptim
    case true
        % Use a log(-x) transform for parameters x = {Qmax_delQ_b, pmax_b, Gmax_b, m_b}
        logNegPars = {'Qmax_delQ_b', 'pmax_b', 'Gmax_b', 'm_b'};
        for i = 1:length(parnames)
            pn = parnames{i};
            if ismember(pn, logNegPars)
                tuneParsTransform.(pn) = @(x) log(-x);
                tuneParsInvTransform.(pn) = @(x) -exp(x);
            else
                tuneParsTransform.(pn) = @(x) x;
                tuneParsInvTransform.(pn) = @(x) x;
            end
        end
        
    case false
        for i = 1:length(parnames)
            pn = parnames{i};
            tuneParsTransform.(pn) = @(x) x;
            tuneParsInvTransform.(pn) = @(x) x;
        end
end

FixedParams.tuneParsTransform = tuneParsTransform;
FixedParams.tuneParsInvTransform = tuneParsInvTransform;


% Assign parameter bounds to workspace -- transformed if required
% switch rescaleForOptim
%     case true
%         for i = 1:length(logNegPars)
%             if ismember(logNegPars{i}, parnames)
%                 func = tuneParsTransform.(logNegPars{i});
%                 ind = strcmp(parnames, logNegPars{i});
%                 % Add small number to any zeros to be log-transformed --
%                 % cannot have infinities.
%                 if ub(ind) == 0 && lb(ind) < 0
%                     ub(ind) = -abs(ub(ind) - lb(ind)) / 1e4;
%                 end
%                 lb(ind) = func(lb(ind));
%                 ub(ind) = func(ub(ind));
%             end
%         end
% end
% assignin('caller', 'npars', npars)
% assignin('caller', 'boundsLower', lb)
% assignin('caller', 'boundsUpper', ub)

switch rescaleForOptim
    case true
        for i = 1:length(logNegPars)
            if ismember(logNegPars{i}, parnames)
                func = tuneParsTransform.(logNegPars{i});
                ind = strcmp(parnames, logNegPars{i});
                bb = [lb(ind) ub(ind)];
                % Add small number to any zeros to be log-transformed --
                % cannot have infinities.
                if bb(2) == 0 && bb(1) < 0
                    bb(2) = -abs(diff(bb)) / 1e4;
                end
                bb = sort(func(bb));
                lb(ind) = bb(1);
                ub(ind) = bb(2);
            end
        end
end
assignin('caller', 'npars', npars)
assignin('caller', 'boundsLower', lb)
assignin('caller', 'boundsUpper', ub)




%% Select cost function.
% There's a few options for the cost function. Not yet sure which is the
% best... Hellinger2_groupWaterOrigin is most defensible as it makes fewest
% assumptions
[~, ~, costFunctionChoices] = costFunction();
% costFunctionChoices = { ...
%     'LSS', ...
%     'RMS', ...
%     'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal', ...
%     'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet', ...
%     'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet', ...
%     'N_LN-Dir_groupWaterOrigin', ...
%     'Hellinger_groupWaterOrigin', ...
%     'Hellinger2_groupWaterOrigin', ...
%     'Hellinger_MVN_groupWaterOrigin'
%     };

if ~exist('costFunctionType', 'var')
    % costFunctionChoices should be given shorter names...
    %     costFunctionType = costFunctionChoices{4}; % select cost function
    %     costFunctionType = costFunctionChoices{6}; % select cost function
    costFunctionType = 'meanCDFdist_Hellinger';
    if ~ismember(costFunctionType, costFunctionChoices)
        costFunctionType = costFunctionChoices{1};
    end
end
FixedParams.costFunction = costFunctionType;
assignin('caller', 'costFunctionLabel', costFunctionType)

%% Select optimising algorithm (so far only ga is available)
optimiserChoices = {'ga','muga'};
if ~exist('optimiser', 'var')
    optimiser = optimiserChoices{1};
end
FixedParams.optimiser = optimiser;
optimise = str2func(optimiser);
assignin('caller', 'optimise', optimise) % assign optimising algorithm to the workspace


%% Parameters of optimising algorithm

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

% Halt integrations at each trajectory's final sampling events to save time
Forc.integrateFullTrajectory = false; % Integrating full trajectories not required for optimisation


%% Filter forcing- and fitting-data

if ~exist('fitTrajectories', 'var')
    fitTrajectories = [];
end

[Forc, Data] = filterInputByOrigin(Forc, Data, 'fitTrajectories', fitTrajectories);

