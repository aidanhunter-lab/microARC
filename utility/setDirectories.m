function x = setDirectories(varargin)
% Output a structure storing all required directories.

% Fields of x may be input as varargin name-value pairs. Although, asides
% from the choice of bioModel, it is probably better to directly modify
% this setDirectories.m script.
extractVarargin(varargin)

% Choose model type
if ~exist('bioModel', 'var')
    bioModelOptions = {'singlePredatorClass', 'multiplePredatorClasses', 'mixotrophy'};
    bioModel = bioModelOptions{2}; % Only singlePredatorClass and multiplePredatorClass available so far
end
x.bioModel = bioModel;

% Forcing data model
if ~exist('forcModel', 'var')
    forcModel = 'sinmod';
end
x.forcModel = forcModel;

% Forcing name (subfolder)
if ~exist('forcName', 'var')
    forcName = 'FramStrait_dep0-100';
end
x.forcName = forcName;

% Experiment name
if ~exist('expName', 'var')
  expName = 'AWI-Hausgarten';  
end
x.expName = expName;

% Forcing data directories
if ~exist('forcDir', 'var')
    forcDir = fullfile('DATA', 'particle_trajectories', forcModel, forcName);
end
x.forcDir = forcDir;
x.forcDummy = 'particles_MODEL_DOMAIN_EXPERIMENT-YEAR_YEAR_t*iSub03.mat';

if ~exist(forcDir, 'dir')
    mkdir(forcDir)
end

% In-situ data directory
if ~exist('obsDir', 'var')
    obsDir = fullfile('DATA', 'AWI_Hausgarten');
end
x.obsDir = obsDir;

% Set/create results directory
if ~exist('results', 'dir')
    mkdir('results')
    addpath results
end
resultsDir = fullfile('results', bioModel);
if ~exist(eval('resultsDir'), 'dir')
    mkdir('results', bioModel)
    addpath(resultsDir)
end
x.resultsDir = resultsDir;

% Output plot directory
plotDir = fullfile(eval('resultsDir'), 'plots');
if ~exist(eval('plotDir'), 'dir')
    mkdir(eval('resultsDir'), 'plots')
    addpath(plotDir)
end
x.plotDir = plotDir;

% Model parameters (default or initial values)
if ~exist('parFile', 'var')
    parFile = 'parameterInitialValues_1.mat';
end

initialParamDir = fullfile(resultsDir, parFile);
initialsExist = ~isempty(parFile) && exist(initialParamDir, 'file');
if ~initialsExist
    initialParamDir = [];
end
% paramName = ['parametersInitialValues_' x.bioModel '.mat'];
x.initialParamDir = initialParamDir;
