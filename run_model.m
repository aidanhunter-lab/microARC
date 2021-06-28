%% Size-structured 1D NPZD model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate model output from default or saved parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Set up model

% Refresh the workspace
clear; clc; close all; delete(gcp('nocreate'))
% Include all subdirectories within search path
addpath(genpath(fileparts(which('run_model'))))

% Store folders/filenames of data and saved parameters
Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', []);
display(Directories)

% Use default set-up to load and organise data and initialise model parameters.
% Set-up may be modified here by passing some extra arguments as name-value
% pairs (preferable), or directly modified within modelSetUp.m
[Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
    'displayAllOutputs', true); % default set-up -- no plots

% Useful name-value pairs for modelSetUp.m include: useTraj, numTraj, ESDmin, ESDmax, nsizes

% modelSetup.m may also produce plots -- set to 'true' any name-value pair as shown below

% [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
%     'plotCellConcSpectra', true, 'plotBioVolSpectra', true, ...
%     'plotSizeClassIntervals', true, ...
%     'trajectoryPlot', true, 'dendrogramPlot', true, ...
%     'plotScalarData', true, 'plotSizeData', true, 'plotAllData', true, ...
%     'displayForc', true, 'displayData', true);

% Parameters may be modified using name-value pairs in updateParameters.m
% Params = updateParameters(Params, FixedParams, 'pmax_a', 20, 'aP', 0.01);

% Input data (forcing trajectories and fitting data) may be filtered by the
% origin of particle trajectories (Atlantic or Arctic)
[Forc, Data] = filterInputByOrigin(Forc, Data, 'fitTrajectories', 'Atlantic');


%% Integrate

% Initialise state variables.
% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
v0 = initialiseVariables(FixedParams, Params, Forc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [size, depth, nutrient]
%                      organic matter      [type, depth, nutrient]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Run the model
tic; disp('.. started at'); disp(datetime('now'))
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ... 
    FixedParams.odeIntegrator, FixedParams.odeOptions);
toc

display(out)
display(auxVars)





