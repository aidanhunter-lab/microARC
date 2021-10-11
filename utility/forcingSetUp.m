function Forc = forcingSetUp(Directories, FixedParams, varargin)

% Derive origin (Arctic or Atlantic) of all particles in Forcing data.
% Save this extra grouping field into the forcing data so that this
% time-consuming function inly needs called once.

extractVarargin(varargin)
extractStruct(Directories)

if ~exist('useTraj', 'var')
    useTraj = [];
end

% Specify years of available forcing data
if ~exist('years', 'var')
    years = 2017:2018;
end

if exist('year', 'var')
    if ~ismember(eval('year'), FixedParams.years)
        error('Chosen year does not match value in struct FixedParams')
    end
else
    year = FixedParams.years;
end

% Load and extract relevant forcing data from physical model 'forcModel'.
F = loadForcing(Directories, years, useTraj);

% Interpolate forcing data over modelled depth layers, and combine multiple
% years of forcing data into a single structure using prepareForcing.m.
% An extra few useful metrics are also calculated here.
Forc = prepareForcing(F, FixedParams);

fileName = ['particleGroupings_' num2str(year) '.mat'];
loadGroups = exist(fileName, 'file');
if loadGroups
    loadedGroups = load(fullfile(forcDir,fileName));
    particleGrouping = loadedGroups.particleGrouping;
end
match = length(Forc.iTraj) == height(particleGrouping); % do loaded values match forcing data?
if loadGroups && match
    Forc.waterMass = particleGrouping.waterMass';
else
    error('The saved table of particle groupings does not match the forcing data!')
end

Forc.nTraj = length(Forc.iTraj);

