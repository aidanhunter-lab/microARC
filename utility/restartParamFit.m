function restartParamFit(fileName, optimOptions, varargin)

% Load output from saved optimisation run and use it to restart the
% optimising algorithm from the previous best position.
% Loaded variables may or may not be assigned to the caller workspace
% depending on varargin options.

if exist(fileName, 'file') ~= 2
    error('Argument "fileName" does not exist as a .mat file.')
end

extractVarargin(varargin)

if ~exist('ni', 'var')
    % New initial values for state variables
    ni = true; % by default, ni = true does not load previous v0 values
end
if ~exist('fp', 'var')
    % Full population of parameters from previous run
    fp = false; % by default, fp = false uses only the best param set found so far, and otherwise randomly initialises parameters 
end

if ~exist('loadMainStructs', 'var')
    % Use loaded values for 4 main data/param structs -- overwrite current workspace?
    loadMainStructs = true;
end
if ~exist('loadData', 'var')
    loadData = true;
end
if ~exist('loadForc', 'var')
    loadForc = true;
end
if ~exist('loadParams', 'var')
    loadParams = true;
end
if ~exist('loadFixedParams', 'var')
    loadFixedParams = true;
end
if ~exist('loadParamBounds', 'var')
    % Use loaded values parameter bounds -- overwrite current workspace?
    loadParamBounds = true;
end

% Load saved outputs
[~, results, ~, ~, boundsLower, boundsUpper, Data, Forc, ... 
    FixedParams, Params, v0] = loadOptimisationRun(fileName);

% Update optimiser initialisation from loaded outputs
populationOld = results.populationHistory(:,:,end);
scoresOld = results.scoreHistory(:,end);
switch fp
    case true
        optimiserOptions = results.optimiserOptions;
        optimiserOptions.InitialPopulationMatrix = populationOld;
        optimiserOptions.InitialScoresMatrix = scoresOld;
    case false
        optimiserOptions = optimOptions;
        optimiserOptions.InitialPopulationMatrix = results.optPar_searchSpace;
end
optimiserOptions.MaxGenerations = optimOptions.MaxGenerations;

% Assign variables to workspace
assignin('caller', 'results', results)
assignin('caller', 'optimiserOptions', optimiserOptions)

switch ni, case false
    assignin('caller', 'v0', v0)
end

switch loadMainStructs, case true
    switch loadData, case true, assignin('caller', 'Data', Data); end
    switch loadForc, case true, assignin('caller', 'Forc', Forc); end
    switch loadParams, case true, assignin('caller', 'Params', Params); end
    switch loadFixedParams, case true, assignin('caller', 'FixedParams', FixedParams); end
end

switch loadParamBounds, case true
    assignin('caller', 'boundsLower', boundsLower)
    assignin('caller', 'boundsUpper', boundsUpper)
end




