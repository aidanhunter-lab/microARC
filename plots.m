%% Size-structured 1D NPZD model

%~~~~~~~~~~~~~~~~~~~
% Plot model outputs
%~~~~~~~~~~~~~~~~~~~

%% Set up model & generate outputs

% Refresh the workspace
clear; close all; delete(gcp('nocreate')); clc
% Include all subdirectories within search path
addpath(genpath(fileparts(which('fit_parameters'))))

% Store folders/filenames of data and saved parameters
Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', []);
% Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
%     'parFile', 'parameterInitialValues_1.mat');
display(Directories)

% Load saved outputs from optimisation runs or choose default set-up
loadFittedParams = true; % use output saved from optimisation run?
fileName = 'fittedParameters';  % saved parameters file name
% tag = '1';                      % and identifying tag
% tag = 'N_LN-Dir_groupWaterOrigin';

% tag = 'Hellinger2_groupWaterOrigin_Atlantic';
% tag = 'Hellinger2_groupWaterOrigin_Arctic';
% tag = 'Hellinger2_groupWaterOrigin_fittedMortalityIntercept_Arctic';
% tag = 'Hellinger3_groupWaterOrigin_Atlantic';
tag = 'IQD_Hellinger_groupWaterOrigin_Atlantic';
tag = 'meanCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj';


% tag = 'Hellinger_MVN_groupWaterOrigin';10
% tag = 'Hellinger_MVN_groupWaterOrigin_2';

fileName = fullfile(Directories.resultsDir, ...
    [fileName '_' tag]);

switch loadFittedParams
    case true
        % Load stored results
        [~, results, parNames, optPar, boundsLower, boundsUpper, Data, Forc, ...
            FixedParams, Params, v0] = loadOptimisationRun(fileName);
        displayFittedParameters(results)
        Params = updateParameters(Params, FixedParams, optPar); % ensure Params struct contains the best-fitting parameter set        
        % Parameters may be modified using name-value pairs in updateParameters.m
        % Params = updateParameters(Params, FixedParams, 'pmax_a', 2, 'aP', 0.01);
        
        % As saved results may be based on trajectories originating fom the
        % Arctic or Atlantic, call modelSetUp to generate the unfiltered
        % forcing and fitting data
        [Forc0, ~, ParamsDefault, Data0] = modelSetUp(Directories, 'numTraj', 1);
        
    case false % Use default model set-up if fitted outputs are not loaded
        [Forc0, FixedParams, Params, Data0] = modelSetUp(Directories);
        
        % modelSetup.m may also produce plots -- set to 'true' any name-value pair as shown below
        
        % [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
        %     'plotCellConcSpectra', true, 'plotBioVolSpectra', true, ...
        %     'plotSizeClassIntervals', true, ...
        %     'trajectoryPlot', true, 'dendrogramPlot', true, ...
        %     'plotScalarData', true, 'plotSizeData', true, 'plotAllData', true, ...
        %     'displayForc', true, 'displayData', true);
end

% Run model over entire trajectories?
Forc.integrateFullTrajectory = true;
% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end
% Initialise variables
if ~exist('v0', 'var') || ~isnumeric(v0)
    % If initial condition v0 has not been loaded then create initials.
    % NOTE: initialiseVariables.m uses the quota parameters to initialise 
    % the plankton state variables. Thus, loaded v0 values may differ from 
    % those here generated using "best-fitting" quota params...
    v0 = initialiseVariables(FixedParams, ParamsDefault, Forc);
end
% v00 = initialiseVariables(FixedParams, Params, Forc0); % create initials for full data (Forc0)
v00 = initialiseVariables(FixedParams, ParamsDefault, Forc0); % create initials for full data (Forc0)


clear out auxVars modData0 modData
% Generate model outputs
% using all trajectories
tic; disp('.. started at'); disp(datetime('now'))
[out0, auxVars0] = integrateTrajectories(FixedParams, Params, Forc0, v00, ... 
    FixedParams.odeIntegrator, FixedParams.odeOptions);
toc
% and only trajectories used to fit the model
tic; disp('.. started at'); disp(datetime('now'))
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ... 
    FixedParams.odeIntegrator, FixedParams.odeOptions);
toc


% Generate modelled equivalents of the data
% modData = matchModOutput2Data(out, auxVars, Data, FixedParams);
modData0 = matchModOutput2Data(out0, auxVars0, Data0, FixedParams);
modData = matchModOutput2Data(out, auxVars, Data, FixedParams);

[cost, costComponents] = costFunction('label', FixedParams.costFunction, ...
    'Data', Data, 'modData', modData);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Plots

% The above outputs are now used as arguments to various plotting 
% functions stored in utility/plottingFunctions/...

folder = Directories.plotDir; % save plots here


%% Map plots

close all

% Particle trajectories from physical model

projection = 'lambert';
% projection = 'miller';
alphaLine = 0.15;

% plt_mapTraj = plot_trajectoryMap(Directories, Forc0, 'projection', projection, ...
%     'alphaLine', alphaLine);
plt_mapTraj = plot_trajectoryMap(Directories, Forc0, 'projection', projection, ...
    'alphaLine', alphaLine, 'Data', Data0);
% plt_mapTraj = plot_trajectoryMap(Directories, Forc0, 'projection', projection, ...
%     'alphaLine', alphaLine, 'Data', Data0, 'labelSampleArea', true);

% plt_mapTraj = plot_trajectoryMap(Directories, Forc0, 'projection', projection, ...
%     'alphaLine', alphaLine, 'Data', Data0, 'flowArrows', true);


plt_mapTraj.Units = 'inches';
plt_mapTraj.Position = [0 0 10 10];


switch save, case true
    if exist('plt_mapTraj', 'var') && isvalid(plt_mapTraj)
        filename = 'map_trajectories.png';
        print(plt_mapTraj, fullfile(folder, filename), '-r300', '-dpng');
    end
end


% In-situ data

projection = 'lambert';
% projection = 'miller';
% projection = 'utm';
alphaPoint = 0.5;
pointSize = 9;
pieSize = 0.02;
lonSpace = 0.4;
latSpace = 0.1;
colourByDataType = true;
showWaterOrigin = true;
polygonAlpha = 0.5;
polygonSmooth = false;
polygonExpand = 0;
omitMapGrid = false; % do not plot map coords -- instead surround data points with polygon used to show sample area in the trajectory map
includeLegend = true;
legendPosition = 'southeast';

plt_mapData = plot_dataMap(Directories, Data0, 'projection', projection, ...
    'alphaPoint', alphaPoint, 'pointSize', pointSize, 'lonSpace', lonSpace, ... 
    'latSpace', latSpace, 'Forc', Forc0, ... 
    'showWaterOrigin', showWaterOrigin, 'polygonAlpha', polygonAlpha, ...
    'polygonExpand', polygonExpand, 'polygonSmooth', polygonSmooth, ...
    'colourByDataType', colourByDataType, 'pieSize', pieSize, ...
    'includeLegend', includeLegend, 'legendPosition', legendPosition, ...
    'omitMapGrid', omitMapGrid);

plt_mapData.Units = 'inches';
plt_mapData.Position = [0 0 6 4];

switch save, case true
    if exist('plt_mapData', 'var') && isvalid(plt_mapData)
        filename = 'map_shipData.png';
        print(plt_mapData, fullfile(folder, filename), '-r300', '-dpng');
    end
end


% Combine the above maps
plt_combineMaps = figure;
plt_combineMaps.Units = 'inches';
plt_combineMaps.Position = [0 0 8 8];

% axes('position', [0.2, 0.2, 0.75, 0.75])
axes('position', [0.15, 0.15, 0.85, 0.85])

projection = 'lambert';
alphaLine = 0.15;
includeLegend = true;
legendPosition = 'east';
legendTextSize = 12;
legendTitle = 'Water origin';
legendTitleSize = 12;
polygonLineWidth = 3;
stripedBorder = false;

plot_trajectoryMap(Directories, Forc0, 'projection', projection, ...
    'alphaLine', alphaLine, 'Data', Data0, 'newPlot', false, ...
    'includeLegend', includeLegend, 'legendPosition', legendPosition, ...
    'legendTitle', legendTitle, 'legendTitleSize', legendTitleSize, ...
    'legendTextSize', legendTextSize, 'polygonLineWidth', polygonLineWidth, ...
    'stripedBorder', stripedBorder);

axes('position', [0.1, 0.12, 0.5, 0.5])

projection = 'lambert';
alphaPoint = 0.5;
pointSize = 9;
pieSize = 0.02;
lonSpace = 0.05;
latSpace = 0.05;
colourByDataType = true;
showWaterOrigin = true;
trimPolygons = true; % shape the water-origin polygons to fit neatly into the full area polygon
polygonAlpha = 1; % polygons cannot be transparent because the underlying map shows though
colSat = 0.6; % reduce colour saturation to emulate the transparent colours
polygonSmooth = false;
polygonExpand = 0;
legendPosition = 'west';
legendTitle = 'Data';
omitMapGrid = true; % do not plot map coords -- instead surround data points with polygon used to show sample area in the trajectory map
fullAreaPolygon = true; % draw polygon matching that used in the trajectory map
polygonLineWidth = 3;

plot_dataMap(Directories, Data0, 'projection', projection, ...
    'alphaPoint', alphaPoint, 'pointSize', pointSize, 'lonSpace', lonSpace, ... 
    'latSpace', latSpace, 'Forc', Forc0, ... 
    'showWaterOrigin', showWaterOrigin, 'polygonAlpha', polygonAlpha, ...
    'polygonExpand', polygonExpand, 'polygonSmooth', polygonSmooth, ...
    'colourByDataType', colourByDataType, 'pieSize', pieSize, ...
    'includeLegend', includeLegend, 'legendPosition', legendPosition, ... 
    'omitMapGrid', omitMapGrid, 'colSat', colSat, ... 
    'fullAreaPolygon', fullAreaPolygon, 'polygonLineWidth', polygonLineWidth, ... 
    'trimPolygons', trimPolygons, 'legendTitle', legendTitle, ... 
    'legendTitleSize', legendTitleSize, 'legendTextSize', legendTextSize, ...
    'stripedBorder', stripedBorder, 'newPlot', false);

switch save, case true
    if exist('plt_combineMaps', 'var') && isvalid(plt_combineMaps)
        filename = 'map_shipDataAndTrajectories.png';
        print(plt_combineMaps, fullfile(folder, filename), '-r300', '-dpng');
    end
end





%% Model fit to data

save = true;

%~~~~~~~~~
% Raw data
%~~~~~~~~~

% Scalar data
plt_scalar = figure;
plt_scalar.Units = 'inches';
plt_scalar.Position = [0 0 8 6];
% axesTextSize = 14;
% legendTextSize = 13;
legendTitle = 'Water origin';
% legendTitleSize = 13;
pointAlpha = 0.4;
subplot(2,2,1)
plot_rawData('scalar', 'N', Data0, 'pointAlpha', pointAlpha, ...
    'includeLegend', true, 'legendTitle', legendTitle);
subplot(2,2,2)
plot_rawData('scalar', 'chl_a', Data0, 'pointAlpha', pointAlpha);
subplot(2,2,3)
plot_rawData('scalar', 'PON', Data0, 'pointAlpha', pointAlpha);
subplot(2,2,4)
plot_rawData('scalar', 'POC', Data0, 'pointAlpha', pointAlpha);


% Size data - spectra
plt_spectra = figure;
plt_spectra.Units = 'inches';
plt_spectra.Position = [0 0 8 6];
subplot(2,1,1)
plot_rawData('sizeSpectra', 'CellConc', Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest');
subplot(2,1,2)
plot_rawData('sizeSpectra', 'BioVol', Data0);

% Size data - spectra multipanel plots
Type = 'CellConc';
plt_spectra2 = figure;
plt_spectra2.Units = 'inches';
plt_spectra2.Position = [0 0 16 12];
subplot(2,2,1)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupAutotroph', true);
subplot(2,2,3)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupHeterotroph', true);
subplot(2,2,2)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupAtlantic', true);
subplot(2,2,4)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupArctic', true);

Type = 'BioVol';
plt_spectra3 = figure;
plt_spectra3.Units = 'inches';
plt_spectra3.Position = [0 0 16 12];
subplot(2,2,1)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupAutotroph', true);
subplot(2,2,3)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupHeterotroph', true);
subplot(2,2,2)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupAtlantic', true);
subplot(2,2,4)
plot_rawData('sizeSpectra', Type, Data0, 'includeLegend', true, ...
    'legendPosition', 'southwest', 'groupArctic', true);

% Size data - binned
plt_binned = figure;
plt_binned.Units = 'inches';
plt_binned.Position = [0 0 8 6];
subplot(2,1,1)
plot_rawData('sizeBinned', 'CellConc', Data0, 'pointAlpha', 0.5, ... 
    'includeLegend', true, 'legendPosition', 'south');
subplot(2,1,2)
plot_rawData('sizeBinned', 'BioVol', Data0, 'pointAlpha', 0.5);

switch save, case true
    % scalar data
    if exist('plt_scalar', 'var') && isvalid(plt_scalar)
        filename = 'data_scalar.png';
        print(plt_scalar, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('plt_spectra', 'var') && isvalid(plt_spectra)
        filename = 'data_sizeSpectra.png';
        print(plt_spectra, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('plt_spectra2', 'var') && isvalid(plt_spectra2)
        filename = 'data_sizeSpectra_groupedCellConc.png';
        print(plt_spectra2, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('plt_spectra3', 'var') && isvalid(plt_spectra3)
        filename = 'data_sizeSpectra_groupedBioVol.png';
        print(plt_spectra3, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('plt_binned', 'var') && isvalid(plt_binned)
        filename = 'data_sizeBinned.png';
        print(plt_binned, fullfile(folder, filename), '-r300', '-dpng');
    end
end

%~~~~~~~~~~~~~~~~~~
% Standardised data
%~~~~~~~~~~~~~~~~~~

close all

plt_stnd = figure;
plt_stnd.Units = 'inches';
plt_stnd.Position = [0 0 16 6];
subplot(2,4,1)
plot_standardisedData('scalar', 'N', Data0, 'covariate', 'Depth', ... 
    'pointAlpha', 0.4, 'densityCurve', true, ...
    'includeLegend', false, 'legendPosition', 'west');
subplot(2,4,5)
plot_standardisedData('scalar', 'N', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot(2,4,2)
plot_standardisedData('scalar', 'PON', Data0, 'covariate', 'Depth', ... 
    'pointAlpha', 0.4, 'densityCurve', true);
subplot(2,4,6)
plot_standardisedData('scalar', 'PON', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot(2,4,3)
plot_standardisedData('scalar', 'POC', Data0, 'covariate', 'Depth', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot(2,4,7)
plot_standardisedData('scalar', 'POC', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot(2,4,4)
plot_standardisedData('scalar', 'chl_a', Data0, 'covariate', 'Depth', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot(2,4,8)
plot_standardisedData('scalar', 'chl_a', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);


switch save, case true
    if exist('plt_stnd', 'var') && isvalid(plt_stnd)
        filename = 'data_scalarStandardised.png';
        print(plt_stnd, fullfile(folder, filename), '-r300', '-dpng');
    end
end

%~~~~~~~~~~~~~~~~~~~
% Model fits to data
%~~~~~~~~~~~~~~~~~~~

close all

% Plot representations of the cost function distance metrics.

% Scalar (nutrient) data
plt_scalarCDFs = figure;
plt_scalarCDFs.Units = 'inches';
plt_scalarCDFs.Position = [0 0 12 8];
markerSize = 5;
showTrueDataCDF = true;
colDat = [0, 0, 0];
colMod = [0, 1, 0];
colMod_ = rgb2hsv(colMod);
colMod_(2) = 0.5;
colMod = hsv2rgb(colMod_);
plotTitle = false;
legendPosition = 'southeast';
itraj = [];
% itraj = 1;

subplot(2,2,1)
varLabel = 'N';
yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
plot_compareCDFs(yobs, ymod, varLabel, 'showTrueDataCDF', showTrueDataCDF, ...
    'plotTitle', plotTitle, 'legendPosition', legendPosition, ...
    'markerSize', markerSize, 'colDat', colDat, 'colMod', colMod, ...
    'itraj', itraj);
subplot(2,2,2)
plotLegend = 'false';
varLabel = 'PON';
yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
plot_compareCDFs(yobs, ymod, varLabel, 'showTrueDataCDF', showTrueDataCDF, ...
    'plotTitle', plotTitle, 'plotLegend', plotLegend, ...
    'markerSize', markerSize, 'colDat', colDat, 'colMod', colMod, ...
    'itraj', itraj);
subplot(2,2,4)
varLabel = 'POC';
yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
plot_compareCDFs(yobs, ymod, varLabel, 'showTrueDataCDF', showTrueDataCDF, ...
    'plotTitle', plotTitle, 'plotLegend', plotLegend, ...
    'markerSize', markerSize, 'colDat', colDat, 'colMod', colMod, ...
    'itraj', itraj);
subplot(2,2,3)
varLabel = 'chl_a';
yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
plot_compareCDFs(yobs, ymod, varLabel, 'showTrueDataCDF', showTrueDataCDF, ...
    'plotTitle', plotTitle, 'plotLegend', plotLegend, ...
    'markerSize', markerSize, 'colDat', colDat, 'colMod', colMod, ...
    'itraj', itraj);
titleText = 'Data distributions vs modelled equivalents';
sgtitle(titleText)


% Size data
plt_sizePMFs = figure;
plt_sizePMFs.Units = 'inches';
plt_sizePMFs.Position = [0 0 12 6];

% itraj = 1; % choose a single trajectory from those used to fit the model, or set itraj = [] to display all model realisations
itraj = [];
barplot = true;
% Vector (size) data
varLabel = 'BioVol';
waterMass = 'Atlantic';
ind0 = strcmp(Data.size.dataBinned.Variable, varLabel) & ... 
    strcmp(Data.size.dataBinned.waterMass, waterMass);
% autotrophs
subplot(1,2,1)
trophicLevel = 'autotroph';
xscale = round(FixedParams.PPdia, 2, 'significant');
ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, trophicLevel);
yobs = Data.size.dataBinned.Value(ind);
ymod = modData.size.Value(ind,:);
plot_comparePMFs(yobs, ymod, varLabel, 'waterMass', waterMass, ... 
    'trophicLevel', trophicLevel, 'itraj', itraj, 'xscale', xscale, ...
    'barplot', barplot);

% heterotrophs
subplot(1,2,2)
trophicLevel = 'heterotroph';
ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, trophicLevel);
yobs = Data.size.dataBinned.Value(ind);
ymod = modData.size.Value(ind,:);
plotLegend = 'false';
plot_comparePMFs(yobs, ymod, varLabel, 'waterMass', waterMass, ... 
    'trophicLevel', trophicLevel, 'itraj', itraj, 'xscale', xscale, ...
    'barplot', barplot, 'plotLegend', plotLegend);


% Summary plots displaying model fit to data
logPlot = true; % for scalar data choose logPlot = true or false
pltChl = plot_fitToData('chl_a', Data, modData, logPlot); pause(0.25)
pltPON = plot_fitToData('PON', Data, modData, logPlot); pause(0.25)
pltPOC = plot_fitToData('POC', Data, modData, logPlot); pause(0.25)
logPlot = false;
pltN = plot_fitToData('N', Data, modData, logPlot); pause(0.25)

logPlot = 'loglog'; % for size spectra data choose logPlot = 'loglog' or 'semilogx'

% Comment out Arctic plots because we're fitting to Atlantic data
pltCellConc_Atl_P = plot_fitToData('CellConc', Data, modData, logPlot, ... 
    'trophicGroup', 'autotroph', 'waterOrigin', 'Atlantic'); pause(0.25)
% pltCellConc_Arc_P = plot_fitToData('CellConc', Data, modData, logPlot, ... 
%     'trophicGroup', 'autotroph', 'waterOrigin', 'Arctic'); pause(0.25)
pltCellConc_Atl_Z = plot_fitToData('CellConc', Data, modData, logPlot, ... 
    'trophicGroup', 'heterotroph', 'waterOrigin', 'Atlantic'); pause(0.25)
% pltCellConc_Arc_Z = plot_fitToData('CellConc', Data, modData, logPlot, ... 
%     'trophicGroup', 'heterotroph', 'waterOrigin', 'Arctic'); pause(0.25)

pltBioVol_Atl_P = plot_fitToData('BioVol', Data, modData, logPlot, ... 
    'trophicGroup', 'autotroph', 'waterOrigin', 'Atlantic'); pause(0.25)
% pltBioVol_Arc_P = plot_fitToData('BioVol', Data, modData, logPlot, ... 
%     'trophicGroup', 'autotroph', 'waterOrigin', 'Arctic'); pause(0.25)
pltBioVol_Atl_Z = plot_fitToData('BioVol', Data, modData, logPlot, ... 
    'trophicGroup', 'heterotroph', 'waterOrigin', 'Atlantic'); pause(0.25)
% pltBioVol_Arc_Z = plot_fitToData('BioVol', Data, modData, logPlot, ... 
%     'trophicGroup', 'heterotroph', 'waterOrigin', 'Arctic'); pause(0.25)



switch save, case true
    % scalar data
    if exist('pltChl', 'var') && isvalid(pltChl)
        filename = 'fitToData_chl.png';
        print(pltChl, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltN', 'var') && isvalid(pltN)
        filename = 'fitToData_DIN.png';
        print(pltN, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltPON', 'var') && isvalid(pltPON)
        filename = 'fitToData_PON.png';
        print(pltPON, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltPOC', 'var') && isvalid(pltPOC)
        filename = 'fitToData_POC.png';
        print(pltPOC, fullfile(folder, filename), '-r300', '-dpng');
    end
    
    % size data
    if exist('pltCellConc_Arc_P', 'var') && isvalid(pltCellConc_Arc_P)
        filename = 'fitToData_CellConcSpectra_Arc_P.png';
        print(pltCellConc_Arc_P, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltCellConc_Arc_Z', 'var') && isvalid(pltCellConc_Arc_Z)
        filename = 'fitToData_CellConcSpectra_Arc_Z.png';
        print(pltCellConc_Arc_Z, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltCellConc_Atl_P', 'var') && isvalid(pltCellConc_Atl_P)
        filename = 'fitToData_CellConcSpectra_Atl_P.png';
        print(pltCellConc_Atl_P, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltCellConc_Atl_Z', 'var') && isvalid(pltCellConc_Atl_Z)
        filename = 'fitToData_CellConcSpectra_Atl_Z.png';
        print(pltCellConc_Atl_Z, fullfile(folder, filename), '-r300', '-dpng');
    end
    
    if exist('pltBioVol_Arc_P', 'var') && isvalid(pltBioVol_Arc_P)
        filename = 'fitToData_BioVolSpectra_Arc_P.png';
        print(pltBioVol_Arc_P, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol_Arc_Z', 'var') && isvalid(pltBioVol_Arc_Z)
        filename = 'fitToData_BioVolSpectra_Arc_Z.png';
        print(pltBioVol_Arc_Z, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol_Atl_P', 'var') && isvalid(pltBioVol_Atl_P)
        filename = 'fitToData_BioVolSpectra_Atl_P.png';
        print(pltBioVol_Atl_P, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol_Atl_Z', 'var') && isvalid(pltBioVol_Atl_Z)
        filename = 'fitToData_BioVolSpectra_Atl_Z.png';
        print(pltBioVol_Atl_Z, fullfile(folder, filename), '-r300', '-dpng');
    end
    
end


% Simpler plots -- boxplots against depth

close all

plt = figure;
subplot(2,2,1)
logPlot = false;
plot_fitToData_depthBoxplot('N', Data, modData, logPlot); pause(0.25)
logPlot = true; % for scalar data choose logPlot = true or false
subplot(2,2,2)
plot_fitToData_depthBoxplot('chl_a', Data, modData, logPlot); pause(0.25)
subplot(2,2,3)
plot_fitToData_depthBoxplot('POC', Data, modData, logPlot); pause(0.25)
subplot(2,2,4)
plot_fitToData_depthBoxplot('PON', Data, modData, logPlot); pause(0.25)
plt.Units = 'inches';
plt.Position = [0 0 8 8];

switch save, case true
    if exist('plt', 'var') && isvalid(plt)
        filename = 'fitToData_depthBoxplot.png';
        print(plt, fullfile(folder, filename), '-r300', '-dpng');
    end
end


%% Fitted parameters

% Display fitted parameters in relation to their bounding values (in the
% table, columns widths shoukd be adjustable).
plt = plot_fittedParameters(results.optPar_summary);

switch save, case true
    if exist('plt', 'var') && isvalid(plt)
        filename = 'fittedParameters.png';
        print(plt, fullfile(folder, filename), '-r300', '-dpng');
    end
end


close all


%% Contour plots -- single trajectories, or grouped by sample event


save = false;

% Choose one or more trajectories -- if multiple are selected then the plot
% will average over them.
sampleEvent = 1;
% All trajectories used for sampleEvent
traj = find(Data.scalar.EventTraj(sampleEvent,:));
% If waterMass is either Atlantic OR Arctic then it may make sense to plot
% average over all trajectories, although there could be unwanted smoothing
% effects...
% Otherwise, if waterMass is a mixture of origins, then group trajectories
% by origin and make separate plots
waterMass = Data.scalar.waterMass{sampleEvent};

switch waterMass
    case {'Atlantic', 'Arctic'}
        plt_Forc = plot_contour_DepthTime('forcing', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_DIN = plot_contour_DepthTime('DIN', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_OM = plot_contour_DepthTime('DOM_POM', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_N = plot_contour_DepthTime('phytoplankton_N', ...
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_Chl = plot_contour_DepthTime('phytoplankton_Chl', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_C = plot_contour_DepthTime('phytoplankton_C', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_N_C = plot_contour_DepthTime('phytoplankton_N_C', ...
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_Chl_N = plot_contour_DepthTime('phytoplankton_Chl_N', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_Z_C = plot_contour_DepthTime('zooplankton_C', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_Z_N = plot_contour_DepthTime('zooplankton_N', ... 
            traj, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
    case 'Arctic/Atlantic'
        if strcmp(waterMass, 'Arctic/Atlantic')
            trajAtlantic = traj(strcmp(Forc.waterMass(traj), 'Atlantic'));
            trajArctic = traj(strcmp(Forc.waterMass(traj), 'Arctic'));
        end
        plt_Forc_Atl = plot_contour_DepthTime('forcing', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ... 
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_Forc_Arc = plot_contour_DepthTime('forcing', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ... 
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_DIN_Atl = plot_contour_DepthTime('DIN', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_DIN_Arc = plot_contour_DepthTime('DIN', ... 
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_OM_Atl = plot_contour_DepthTime('DOM_POM', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_OM_Arc = plot_contour_DepthTime('DOM_POM', ... 
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_N_Atl = plot_contour_DepthTime('phytoplankton_N', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_N_Arc = plot_contour_DepthTime('phytoplankton_N', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_Chl_Atl = plot_contour_DepthTime('phytoplankton_Chl', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_Chl_Arc = plot_contour_DepthTime('phytoplankton_Chl', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_C_Atl = plot_contour_DepthTime('phytoplankton_C', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_C_Arc = plot_contour_DepthTime('phytoplankton_C', ... 
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_N_C_Atl = plot_contour_DepthTime('phytoplankton_N_C', ... 
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_N_C_Arc = plot_contour_DepthTime('phytoplankton_N_C', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_P_Chl_N_Atl = plot_contour_DepthTime('phytoplankton_Chl_N', ...
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_P_Chl_N_Arc = plot_contour_DepthTime('phytoplankton_Chl_N', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_Z_C_Atl = plot_contour_DepthTime('zooplankton_C', ...
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_Z_C_Arc = plot_contour_DepthTime('zooplankton_C', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
        plt_Z_N_Atl = plot_contour_DepthTime('zooplankton_N', ...
            trajAtlantic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Atlantic');
        plt_Z_N_Arc = plot_contour_DepthTime('zooplankton_N', ...
            trajArctic, out, auxVars, FixedParams, Forc, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', 'Arctic');
end

switch save, case true
    
    switch waterMass
        
        case {'Atlantic', 'Arctic'}
        
            if exist('plt_Forc', 'var') && isvalid(plt_Forc)
                filename = ['forcing_data_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_Forc, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_DIN', 'var') && isvalid(plt_DIN)
                filename = ['DIN_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_DIN, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_OM', 'var') && isvalid(plt_OM)
                filename = ['OM_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_OM, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_N', 'var') && isvalid(plt_P_N)
                filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_N, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_Chl', 'var') && isvalid(plt_P_Chl)
                filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_Chl, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_C', 'var') && isvalid(plt_P_C)
                filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_C, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_N_C', 'var') && isvalid(plt_P_N_C)
                filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_N_C, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_P_Chl_N', 'var') && isvalid(plt_P_Chl)
                filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_P_Chl_N, fullfile(folder, filename), '-r300', '-dpng');
            end            
            if exist('plt_Z_C', 'var') && isvalid(plt_Z_C)
                filename = ['zooplankton_C_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_Z_C, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Z_N', 'var') && isvalid(plt_Z_N)
                filename = ['zooplankton_N_sampleEvent' num2str(sampleEvent) '.png'];
                print(plt_Z_N, fullfile(folder, filename), '-r300', '-dpng');
            end
            
        case  'Arctic/Atlantic'
            
            if exist('plt_Forc_Atl', 'var') && isvalid(plt_Forc_Atl)
                filename = ['forcing_data_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_Forc_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Forc_Arc', 'var') && isvalid(plt_Forc_Arc)
                filename = ['forcing_data_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_Forc_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end

            if exist('plt_DIN_Atl', 'var') && isvalid(plt_DIN_Atl)
                filename = ['DIN_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_DIN_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_DIN_Arc', 'var') && isvalid(plt_DIN_Arc)
                filename = ['DIN_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_DIN_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end

            if exist('plt_OM_Atl', 'var') && isvalid(plt_OM_Atl)
                filename = ['OM_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_OM_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_OM_Arc', 'var') && isvalid(plt_OM_Arc)
                filename = ['OM_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_OM_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_N_Atl', 'var') && isvalid(plt_P_N_Atl)
                filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_N_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_N_Arc', 'var') && isvalid(plt_P_N_Arc)
                filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_N_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_Chl_Atl', 'var') && isvalid(plt_P_Chl_Atl)
                filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_Chl_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_Chl_Arc', 'var') && isvalid(plt_P_Chl_Arc)
                filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_Chl_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_C_Atl', 'var') && isvalid(plt_P_C_Atl)
                filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_C_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_C_Arc', 'var') && isvalid(plt_P_C_Arc)
                filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_C_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_N_C_Atl', 'var') && isvalid(plt_P_N_C_Atl)
                filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_N_C_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_N_C_Arc', 'var') && isvalid(plt_P_N_C_Arc)
                filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_N_C_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
            if exist('plt_P_Chl_N_Atl', 'var') && isvalid(plt_P_Chl_N_Atl)
                filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_P_Chl_N_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_P_Chl_N_Arc', 'var') && isvalid(plt_P_Chl_N_Arc)
                filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_P_Chl_N_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Z_C_Atl', 'var') && isvalid(plt_Z_C_Atl)
                filename = ['zooplankton_C_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_Z_C_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Z_C_Arc', 'var') && isvalid(plt_Z_C_Arc)
                filename = ['zooplankton_C_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_Z_C_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Z_N_Atl', 'var') && isvalid(plt_Z_N_Atl)
                filename = ['zooplankton_N_sampleEvent' num2str(sampleEvent) '_AtlanticOrigin.png'];
                print(plt_Z_N_Atl, fullfile(folder, filename), '-r300', '-dpng');
            end
            if exist('plt_Z_N_Arc', 'var') && isvalid(plt_Z_N_Arc)
                filename = ['zooplankton_N_sampleEvent' num2str(sampleEvent) '_ArcticOrigin.png'];
                print(plt_Z_N_Arc, fullfile(folder, filename), '-r300', '-dpng');
            end
            
    end
end


%% Time evolution

% These polygon plots should be extended to show water masses of Atlantic
% and of Arctic orgin because some sampling events use trajectories
% originating from both oceans... red for Atlantic blue for Arctic

close all

% Choose event
sampleEvent = 1;
if ~ismember(sampleEvent, 1:Data.scalar.nEvents), warning(['Choose event number within range (1, ' num2str(Data.scalar.nEvents) ')']); end
% trajectory indices
traj = find(Data.scalar.EventTraj(sampleEvent,:));

highlightColour = [1 0 1];
plotOptions = {'forcing', 'DIN', 'organicN', 'organicC', 'phytoplankton_C', ...
    'phytoplanktonStacked', 'phytoZooPlanktonStacked'};

for varIndex = 1:length(plotOptions)
    
    plotVariables = plotOptions{varIndex};
    
    switch plotVariables
        
        case 'forcing'
            % Forcing data
            plt_forcing = figure;
            plt_forcing.Units = 'inches';
            plt_forcing.Position = [0 0 12 12];
            
            % temperature
            subplot(3,1,1)
            plot_timeSeries_trajectoryPolygon('temperature', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            % diffusivity
            subplot(3,1,2)
            plot_timeSeries_trajectoryPolygon('diffusivity', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            % PAR
            subplot(3,1,3)
            plot_timeSeries_trajectoryPolygon('PAR', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'DIN'
            % DIN
            plt_DIN = figure;            
            plt_DIN.Units = 'inches';
            plt_DIN.Position = [0 0 12 8];
            
            subplot(2,1,1)
            plot_timeSeries_trajectoryPolygon('DIN', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,1,2)
            plot_timeSeries_trajectoryPolygon('DIN', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'organicN'
            plt_ON = figure;
            plt_ON.Units = 'inches';
            plt_ON.Position = [0 0 24 8];
            
            subplot(2,2,1)
            plot_timeSeries_trajectoryPolygon('DON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,2)
            plot_timeSeries_trajectoryPolygon('PON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,3)
            plot_timeSeries_trajectoryPolygon('DON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,4)
            plot_timeSeries_trajectoryPolygon('PON', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'organicC'
            plt_OC = figure;
            plt_OC.Units = 'inches';
            plt_OC.Position = [0 0 24 8];

            subplot(2,2,1)
            plot_timeSeries_trajectoryPolygon('DOC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,2)
            plot_timeSeries_trajectoryPolygon('POC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,3)
            plot_timeSeries_trajectoryPolygon('DOC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,4)
            plot_timeSeries_trajectoryPolygon('POC', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'phytoplankton_C'
            plt_P_C = plot_timeSeries_trajectoryPolygon('phytoplankton_C', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'max', 'highlightColour', highlightColour, 'fixedYaxis', false);
            plt_P_C_fixedScale = plot_timeSeries_trajectoryPolygon('phytoplankton_C', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data, ...
                'depth', 'max', 'highlightColour', highlightColour, 'fixedYaxis', true);
            
        case 'phytoplanktonStacked'
            plt_P_C_stacked = plot_timeSeries_trajectoryPolygon('phytoplanktonStacked', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data);

        case 'phytoZooPlanktonStacked'
            plt_P_C_stacked = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
                sampleEvent, traj, out, auxVars, FixedParams, Forc, Data);
    end
end

switch save
    
    case true
        
        if exist('plt_forcing', 'var') && isvalid(plt_forcing)
            filename = ['forcing_data_timeSeriesPolygon_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_forcing, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_DIN', 'var') && isvalid(plt_DIN)
            filename = ['DIN_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_DIN, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_OM', 'var') && isvalid(plt_OM)
            filename = ['OM_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_OM, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_N', 'var') && isvalid(plt_P_N)
            filename = ['phytoplankton_N_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_N, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_Chl', 'var') && isvalid(plt_P_Chl)
            filename = ['phytoplankton_Chl_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_Chl, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_C', 'var') && isvalid(plt_P_C)
            filename = ['phytoplankton_C_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_C, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_N_C', 'var') && isvalid(plt_P_N_C)
            filename = ['phytoplankton_N_C_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_N_C, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_P_Chl_N', 'var') && isvalid(plt_P_Chl_N)
            filename = ['phytoplankton_Chl_N_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_P_Chl_N, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_Z', 'var') && isvalid(plt_Z)
            filename = ['zooplankton_sampleEvent' num2str(sampleEvent) '.png'];
            print(plt_Z, fullfile(folder, filename), '-r300', '-dpng');
        end
end


% Group trajectories of Arctic and of Atlantic origin

save = false;

plt_Atlantic = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
    [], [], out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', 'Atlantic');
plt_Arctic = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
    [], [], out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', 'Arctic');



% axesTextSize = 16;
% titleTextSize = 16;
% legendTextSize = 16;
% legendPointSize = 60;
% plt_Atlantic = plot_timeSeries_trajectoryPolygon2('phytoZooPlanktonStacked', ...
%     [], [], out, auxVars, FixedParams, Forc0, Data0, ...
%     'waterMass', 'Atlantic', 'axesTextSize', axesTextSize, ... 
%     'titleTextSize', titleTextSize, 'legendTextSize', legendTextSize, ...
%     'legendPointSize', legendPointSize);
% plt_Atlantic.Units = 'inches';
% plt_Atlantic.Position = [0 0 10 4];


switch save
    case true
        if exist('plt_Atlantic', 'var') && isvalid(plt_Atlantic)
            filename = 'stackedPolygons_planktonAtlantic.png';
            print(plt_Atlantic, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_Arctic', 'var') && isvalid(plt_Arctic)
            filename = 'stackedPolygons_planktonArctic.png';
            print(plt_Arctic, fullfile(folder, filename), '-r300', '-dpng');
        end
end



% plt = plot_timeSeries_barplot('phytoZooPlankton',sampleEvent,traj,out,FixedParams,Forc,Data);
% 
% switch save
%     case true
%         if exist('plt', 'var') && isvalid(plt)
%             filename = ['barplot_timeSeries_plankton_sampleEvent' num2str(sampleEvent) '.png'];
%             print(plt, fullfile(folder, filename), '-r300', '-dpng');
%         end
% end
% 



% Make animated plots showing fluctuations in abundance, and possibly
% predator cycles, and probably increased variability within the smallest
% size classes

work in progress...

xp = auxVars0.biovolume(1:9,:,:,1);
xp = auxVars0.cellDensity(1:9,:,:,1);
xp = squeeze(sum(xp, 2));
xz = auxVars0.biovolume(10:18,:,:,1);
xz = auxVars0.cellDensity(10:18,:,:,1);
xz = squeeze(sum(xz, 2));

figure
semilogy(1:9,xp(:,180))
hold on
semilogy(1:9,xz(:,180))

hold off





%% Network plots -- fluxes, production

close all

% Feeding fluxes

plt_feedFlux_N = figure;
plt_feedFlux_N.Units = 'inches';
plt_feedFlux_N.Position = [0 0 10 7.5];
subplot(2,1,1)
plot_network('feedingFluxes', 'nitrogen', auxVars0, FixedParams, Forc0, 'Arctic');
subplot(2,1,2)
plot_network('feedingFluxes', 'nitrogen', auxVars0, FixedParams, Forc0, 'Atlantic');



plt_feedFlux_C = figure;
plt_feedFlux_C.Units = 'inches';
plt_feedFlux_C.Position = [0 0 10 7.5];
subplot(2,1,1)
plot_network('feedingFluxes', 'carbon', auxVars0, FixedParams, Forc0, 'Arctic');
subplot(2,1,2)
plot_network('feedingFluxes', 'carbon', auxVars0, FixedParams, Forc0, 'Atlantic');


% Organic matter
plt_OM = figure;
plt_OM.Units = 'inches';
plt_OM.Position = [0 0 10 10];
subplot(2,2,1)
plot_network('OMfluxes', 'carbon', auxVars0, FixedParams, Forc0, 'Arctic');
subplot(2,2,2)
plot_network('OMfluxes', 'nitrogen', auxVars0, FixedParams, Forc0, 'Arctic');
subplot(2,2,3)
plot_network('OMfluxes', 'carbon', auxVars0, FixedParams, Forc0, 'Atlantic');
subplot(2,2,4)
plot_network('OMfluxes', 'nitrogen', auxVars0, FixedParams, Forc0, 'Atlantic');



switch save, case true
        if exist('plt_feedFlux_C', 'var') && isvalid(plt_feedFlux_C)
            filename = 'networkPlot_feedingFlux_carbon.png';
            print(plt_feedFlux_C, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_feedFlux_N', 'var') && isvalid(plt_feedFlux_N)
            filename = 'networkPlot_feedingFlux_nitrogen.png';
            print(plt_feedFlux_N, fullfile(folder, filename), '-r300', '-dpng');
        end
        if exist('plt_OM', 'var') && isvalid(plt_OM)
            filename = 'networkPlot_OMFluxes.png';
            print(plt_OM, fullfile(folder, filename), '-r300', '-dpng');
        end
end




%% Parameter plots -- correlations

parNames = results.parNames;
popHist = results.populationHistory;
costHist = results.scoreHistory;

size(popHist)
size(costHist)

% index parameter sets within 15% of best cost
optCost = min(costHist(:));
costInd = costHist <= optCost + 0.5 * optCost;
costInd = reshape(costInd, [100, 1, 40]);
% costInd = repmat(reshape(costInd, [100, 1, 40]), [1, length(parNames), 1]);


Gmax_a_i = strcmp(parNames, 'Gmax_a');
k_G_i = strcmp(parNames, 'k_G');
m_a_i = strcmp(parNames, 'm_a');

Gmax_a = popHist(:,Gmax_a_i,:);
k_G = popHist(:,k_G_i,:);
m_a = popHist(:,m_a_i,:);

figure
scatter(Gmax_a(costInd), k_G(costInd))
scatter(m_a(costInd), Gmax_a(costInd))


figure
scatter(Gmax_a(:), costHist(:))

figure
Gmax_b_i = strcmp(parNames, 'Gmax_b');
Gmax_b = -exp(popHist(:,Gmax_b_i,:));
scatter(Gmax_b(:), costHist(:))


figure
m_b_i = strcmp(parNames, 'm_b');
m_b = -exp(popHist(:,m_b_i,:));
scatter(m_b(:), costHist(:))




%% Plot groups of trajectories corresponding to each event

% These plots need improved... better to run model over ALL available
% trajectories (not just those used for optimisation), then display output
% given that greater spatial coverage, and optionally show sampling events
% Also, move the plot functions to their own separate script to get rid of
% the outputPlot.m that aggregates them all...

% Choose event
ie = 7;
if ~ismember(ie, 1:Data.scalar.nEvents), warning(['Choose event number within range (1, ' num2str(Data.scalar.nEvents) ')']); end
% trajectory indices
kk = find(Data.scalar.EventTraj(ie,:));

outputPlot('trajectoryLine_LatLong','direction',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','forcing',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','inorganicNutrient',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','DOM_POM',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
outputPlot('trajectoryLine_LatLong','phytoplankton_N',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);
% outputPlot('trajectoryLine_LatLong','zooplankton',ie,kk,out,FixedParams,Forc,auxVars,Data,0.1);


close all

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% %### in progress...
% 
% ie = 1; % sampling event
% kk = find(Data.scalar.EventTraj(ie,:)); % all trajectories for selected event
% k = kk(7); % single trajectory
% 
% % x = out.N(:,:,:,:);
% x = out.N(:,:,:,k);
% 
% xmin = squeeze(min(x));
% xmax = squeeze(max(x));
% % tt = yearday(Forc.t(:,:));
% tt = yearday(Forc.t(:,k));
% n = size(tt, 1);
% 
% lat = Forc.y(:,k); % latitude
% % lat = Forc.y(:,:); % latitude
% % plot colour changing with latitude along line..?
% 
% xmin = reshape(xmin, [1 n]);
% xmax = reshape(xmax, [1 n]);
% tt = reshape(tt, [1 n]);
% lat = reshape(lat, [1 n]);
% z = zeros(size(xmin));
% 
% figure
% surface([tt;tt], [xmin;xmin], [z;z], [lat;lat], ...
%     'facecol', 'no', 'edgecol', 'interp', 'linew', 2)
% hold on
% surface([tt;tt], [xmax;xmax], [z;z], [lat;lat], ...
%     'facecol', 'no', 'edgecol', 'interp', 'linew', 2)
% hold off
% cb = colorbar();
% deg = char(176);
% cb.Label.String = ['latitude (' deg 'N)'];
% xlabel('day of year')
% ylabel(['DIN (mmol N m^{-3})'])
% title('max and min DIN over time')
% 
% 
% colormap(plasma)
% 
% % ###
% % Aggregate trajectories originating from Atlantic and Arctic waters then
% % make summmary plots displaying typical (across-trajectories) results...
% 
% pltN = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, 'DIN', 'minmax');
% 
% pltOM_N = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, 'OM_N', 'max');
% pltOM_C = TimeEvolution_ColByLat_GroupTrajectories(out, Forc, FixedParams, 'OM_C', 'max');
% 
% now do the zooplanton... save phyto till last...
% 

