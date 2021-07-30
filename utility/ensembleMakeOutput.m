function [modData, out, auxVars] = ensembleMakeOutput(pars, ... 
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, varargin)
% Wrapper function for integrateTrajectories.m and matchModOutput2Data.m.
% To reduce memory requirements of outputs, include optional arguments to
% filter the 'auxVars' and 'out' outputs... still to do...

extractVarargin(varargin)

returnExtra_ = {'cellDensity', 'biovolume'}; % extra output needed for the cost function
if ~exist('returnExtra', 'var')
    % more outputs may be returned as extra outputs (not recommended while
    % optimising parameters)
    returnExtra = returnExtra_;
else
    returnExtra = unique([eval('returnExtra'), returnExtra_], 'stable');
end
    
%% Run model

% Set parameter values
Params = updateParameters(Params, FixedParams, pars);

% Integrate
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ...
    odeIntegrator, odeOptions, 'returnExtra', returnExtra);



%% Match model outputs to data

% Extract model output from times and depths matching the data, then 
% transform it using the same functions that standardised the data. The
% resulting 'modData' should have the same structure as, and be directly 
% comparable to, the observations in 'Data'.

modData = matchModOutput2Data(out, auxVars, Data, FixedParams, ...
    'fitToFullSizeSpectra', FixedParams.fitToFullSizeSpectra);

%%

if ~exist('modData', 'var')
    modData = nan;
end
if ~exist('out', 'var')
    out = nan;
end
if ~exist('auxVars', 'var')
    auxVars = nan;
end
