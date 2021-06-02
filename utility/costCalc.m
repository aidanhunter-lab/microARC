function [cost, costComponents, modData, out, auxVars] = costCalc(pars, ... 
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, varargin)
% Returns a scalar 'cost' quantifying model misfit to data, given parameters
% 'pars'.
% May optionally return a struct 'costComponents' containing the cost 
% ascribed to each separate data type; all model outputs, 'out' and 
% 'auxVars'; and the modelled equivalents of the data, 'modData'.

extractVarargin(varargin)

returnExtra_ = {'cellDensity', 'biovolume'}; % extra output needed for the cost function
if ~exist('returnExtra', 'var')
    % more outputs may be returned as extra outputs (not recommended while
    % optimising parameters)
    returnExtra = returnExtra_;
else
    returnExtra = unique([eval('returnExtra'), returnExtra_], 'stable');
end
    
% Cost function type should be given by name-value pair in varargin,
% otherwise the default is used
defaultCostType = 'LSS';
if ~exist('selectFunction', 'var')
    selectFunction = defaultCostType;
end

% Vars = Data.scalar.obsInCostFunction; % Data types used to fit the model
% VarsSize = Data.size.obsInCostFunction;


%% Run model

% Set parameter values
Params = updateParameters(Params, FixedParams, pars);

% % Some initial state variables are specified as functions of parameter
% % values => should be recalculated for each parameter set. If v0 were totally
% % independent of parameter values then it would be absolutely fixed and, to 
% % avoid wasted computation, should be passed as an argument to costFunction.m.
% v0 = initialiseVariables(FixedParams, Params, Forc);

% Integrate
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ...
    odeIntegrator, odeOptions, 'returnExtra', returnExtra);


%% Match model outputs to data

% Extract model output from times and depths matching the data, then 
% transform it using the same functions that standardised the data. The
% resulting 'modData' should have the same structure as, and be directly 
% comparable to, the observations in 'Data'.
modData = matchModOutput2Data(out, auxVars, Data, FixedParams);


%% Cost function

[cost, costComponents] = costFunction('label', selectFunction, ...
    'Data', Data, 'modData', modData);


%%

if ~exist('cost', 'var')
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if ~exist('costComponents', 'var')
    costComponents = nan;
end
if ~exist('modData', 'var')
    modData = nan;
end
if ~exist('out', 'var')
    out = nan;
end
if ~exist('auxVars', 'var')
    auxVars = nan;
end

