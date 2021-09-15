%% Size-structured 1D NPZD model

%~~~~~~~~~~~~~~~~~~~
% Plot model outputs
%~~~~~~~~~~~~~~~~~~~

%% Set up model & generate outputs

% Refresh the workspace
clear; close all; delete(gcp('nocreate')); clc
% Include all subdirectories within search path
addpath(genpath(fileparts(which('plots'))))

% Store folders/filenames of data and saved parameters
% file_params = [];
fileName_params = 'parameterInitialValues';
fileType_params = '.mat';
tag = '_RMS_Hellinger2_Atlantic_aG_sigG_upweightAbnTot';
file_params = fullfile([fileName_params tag fileType_params]);

if ~exist(file_params, 'file'), warning('file_params does not exist -- check file names in results directory.'), file_params = []; end

Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', file_params);
display(Directories)

% Load saved outputs from optimisation runs or choose default set-up
loadFittedParams = true; % use output saved from optimisation run?
% Saved parameters file name and identifying tag
fileName_results = 'fittedParameters';
fileType_results = '.mat';
tag = '_RMS_Hellinger2_Atlantic_aG_sigG_upweightAbnTot';
file_results = fullfile(Directories.resultsDir, ...
    [fileName_results tag fileType_results]);

switch loadFittedParams
    case true
        % Load stored results
        [~, results, parNames, optPar, boundsLower, boundsUpper, Data, Forc, ...
            FixedParams, Params, v0] = loadOptimisationRun(file_results);
        Params = updateParameters(Params, FixedParams, optPar); % update Params struct to store the best-fitting parameter set        
        % Loaded parameters and associated data may be based on (filtered)
        % particle trajectories originating fom the Arctic or Atlantic.
        % Some plots will show all data so call modelSetUp to generate the
        % complete forcing & fitting data sets.
        numTraj = 1; % Single trajectory per sampling event
        chlSampleDepthLimit = inf; % Include all chlorophyll samples -- do not omit the deep samples
        [Forc0, ~, ~, Data0] = modelSetUp(Directories, ... 
            'numTraj', numTraj, 'chlSampleDepthLimit', chlSampleDepthLimit);
    
    case false % Use default model set-up if fitted outputs are not loaded
        numTraj = 1;
        chlSampleDepthLimit = inf;
        [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
            'numTraj', numTraj, 'chlSampleDepthLimit', chlSampleDepthLimit);
        Forc0 = Forc; Data0 = Data;
end

% Run model over entire trajectories?
Forc.integrateFullTrajectory = true;
Forc0.integrateFullTrajectory = Forc.integrateFullTrajectory;

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Initialise variables
if ~exist('v0', 'var') || ~isnumeric(v0)
    % If initial condition v0 has not been loaded then create initials.
    v0 = initialiseVariables(FixedParams, Params, Forc);
end
v00 = initialiseVariables(FixedParams, Params, Forc0); % create initials for full data (Forc0)

% Generate model outputs over trajectories linked to all data (Atlantic & Arctic)
clear out out0 auxVars auxVars0 modData0 modData
tic; disp('.. started at'); disp(datetime('now'))
[out0, auxVars0] = integrateTrajectories(FixedParams, Params, Forc0, v00, ... 
    FixedParams.odeIntegrator, FixedParams.odeOptions);
toc
% and only over the trajectories used to fit the model
tic; disp('.. started at'); disp(datetime('now'))
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ... 
    FixedParams.odeIntegrator, FixedParams.odeOptions);
toc

% Generate modelled equivalents of the data
if ~isfield(FixedParams, 'fitToFullSizeSpectra'), FixedParams.fitToFullSizeSpectra = false; end
modData0 = matchModOutput2Data(out0, auxVars0, Data0, FixedParams, ...
    'fitToFullSizeSpectra', FixedParams.fitToFullSizeSpectra);
modData = matchModOutput2Data(out, auxVars, Data, FixedParams, ...
    'fitToFullSizeSpectra', FixedParams.fitToFullSizeSpectra);

[cost, costComponents] = costFunction('label', FixedParams.costFunction, ...
    'Data', Data, 'modData', modData);
[cost0, costComponents0] = costFunction('label', FixedParams.costFunction, ...
    'Data', Data0, 'modData', modData0);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Plots

% Use the above outputs as arguments to plotting functions stored in utility/plottingFunctions/...

folder = Directories.plotDir; % save plots here

save = true;

%% Display data

% Nutrient & organic matter at depth
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


% Standardised data -- sized for paper
plt_stnd = figure;
plt_stnd.Units = 'inches';
plt_stnd.Position = [0 0 6 9];
nrows = 4; % number of rows excluding legend
ncols = 2;
legh = (1 / 5) * (1 / nrows); % legend height
legFontSize = 11;
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height (all panels and legend)
Pw = 1 - 2 * pex; % total plot width
yex = 0.3; % proportion of panel size devoted to axis labels -- y dimension
xex = 0.2;
pht = (Ph - legh) / nrows; % panel height total
ph = pht - yex * pht; % panel height excluding labels
pwt = Pw / ncols;
pw = pwt - xex * pwt;

subplot('Position', [pex + (pwt - pw), 1 - pex - ph, pw, ph])
plot_standardisedData('scalar', 'N', Data0, 'covariate', 'Depth', ... 
    'pointAlpha', 0.4, 'densityCurve', true, ...
    'includeLegend', false, 'legendPosition', 'west');
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - ph, pw, ph])
plot_standardisedData('scalar', 'N', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + (pwt - pw), 1 - pex - pht - ph, pw, ph])
plot_standardisedData('scalar', 'PON', Data0, 'covariate', 'Depth', ... 
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - pht - ph, pw, ph])
plot_standardisedData('scalar', 'PON', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + (pwt - pw), 1 - pex - 2 * pht - ph, pw, ph])
plot_standardisedData('scalar', 'POC', Data0, 'covariate', 'Depth', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - 2 * pht - ph, pw, ph])
plot_standardisedData('scalar', 'POC', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + (pwt - pw), 1 - pex - 3 * pht - ph, pw, ph]);
plot_standardisedData('scalar', 'chl_a', Data0, 'covariate', 'Depth', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - 3 * pht - ph, pw, ph]);
plot_standardisedData('scalar', 'chl_a', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
% Include legend at bottom
if ismember('Arctic/Atlantic', Data0.scalar.waterMass)
    waterMasses = {'Arctic', 'Arctic & Atlantic', 'Atlantic'};
    ncol = 3;
else
    waterMasses = {'Arctic', 'Atlantic'};    
    ncol = 2;
end
leg = legend([waterMasses, 'data distribution', 'standard normal'], ...
    'Location', 'bestoutside', 'FontSize', legFontSize, ...
    'Orientation', 'horizontal', 'NumColumns', ncol, ...
    'box', 'on');
set(leg, 'Position', [0.25, pex, 0.5, legh])


% Size spectra -- averaged over samples
plt_spectra = figure;
plt_spectra.Units = 'inches';
plt_spectra.Position = [0 0 8 6];
% meanType = 'arithmetic';
meanType = 'geometric';
legendPosition = 'south';
subplot(2,1,1)
plot_rawData('sizeSpectra', 'CellConc', Data0, 'meanType', meanType);%, ... 
%     'includeLegend', true, 'legendPosition', legendPosition);
subplot(2,1,2)
plot_rawData('sizeSpectra', 'BioVol', Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition);


% Size spectra -- multipanel plots
Type = 'CellConc';
plt_spectra2 = figure;
plt_spectra2.Units = 'inches';
plt_spectra2.Position = [0 0 16 12];
legendPosition = 'south';
subplot(2,2,1)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupAutotroph', true);
subplot(2,2,3)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupHeterotroph', true);
subplot(2,2,2)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ...
    'groupAtlantic', true);
subplot(2,2,4)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupArctic', true);

Type = 'BioVol';
plt_spectra3 = figure;
plt_spectra3.Units = 'inches';
plt_spectra3.Position = [0 0 16 12];
subplot(2,2,1)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupAutotroph', true);
subplot(2,2,3)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupHeterotroph', true);
subplot(2,2,2)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupAtlantic', true);
subplot(2,2,4)
plot_rawData('sizeSpectra', Type, Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition, ... 
    'groupArctic', true);


% Integrated size spectra -- binned size data
plt_binned = figure;
plt_binned.Units = 'inches';
plt_binned.Position = [0 0 8 6];
subplot(2,1,1)
plot_rawData('sizeBinned', 'CellConc', Data0, 'pointAlpha', 0.5, ... 
    'includeLegend', true, 'legendPosition', 'south');
subplot(2,1,2)
plot_rawData('sizeBinned', 'BioVol', Data0, 'pointAlpha', 0.5);


switch save, case true
    if exist('plt_scalar', 'var') && isvalid(plt_scalar)
        filename = 'data_scalar.png';
        print(plt_scalar, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('plt_stnd', 'var') && isvalid(plt_stnd)
        filename = 'data_scalarStandardised.png';
        print(plt_stnd, fullfile(folder, filename), '-r300', '-dpng');
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


close all


%% Model fit to data

% Nutrient (scalar) data. Multipanel summary plots for an overview of model
% fit. Include 3 rows to display: (1) ungrouped data; (2) grouped by
% depth; (3) grouped by event.
standard = [true, false];
pn = {'standardised', 'raw'};
Vars = {'N','PON','POC','chl_a'};
nrows = 3; % 1 row per plot type
ncols = length(Vars); % 1 column per data type
colDat = [0, 0, 0];
colMod = [0, 1, 0];
connectors = true;
for i = 1:2
    standardised = standard(i);
    plotName = ['pltFit2Data_' pn{i}];
    assignin('base', plotName, figure)
    set(eval(plotName), 'Units', 'inches')
    set(eval(plotName), 'Position', [0 0 16 12])
    
    % 1st row: ungrouped data
    for ii = 1:ncols
        subplot(nrows, ncols, ii)
        xvar = Vars{ii};
        plot_fitToNutrient_sorted(xvar, Data, modData, ...
            'colDat', colDat, 'colMod', colMod, ...
            'standardised', standardised, 'connectors', connectors);
    end
    
    % 2nd row: grouped by depth -- boxplot
    for ii = 1:ncols
        subplot(nrows, ncols, ii + ncols)
        xvar = Vars{ii};
        plot_fitToNutrient_depth(xvar, Data, modData, ...
            'colDat', colDat, 'colMod', colMod, ...
            'standardised', standardised);
    end
    
    % 3rd row: grouped by sample event -- boxplot
    for ii = 1:ncols
        subplot(nrows, ncols, ii + 2 * ncols)
        xvar = Vars{ii};
        plot_fitToNutrient_event(xvar, Data, modData, ...
            'colDat', colDat, 'colMod', colMod, ...
            'standardised', standardised);
    end
    sgtitle(['Model fit to ' pn{i} ' data'])
end

% Grouped by data type
standardised = true;
ap = '';
switch standardised, case true, ap = '_standardised'; end
nrows = 1;
ncols = 3;
for i = 1:length(Vars)
    xvar = Vars{i};
    xLab = xvar;
    switch xvar
        case 'N', xLab = 'DIN';
        case 'chl_a', xLab = 'chl a';
    end    
    plotName = ['pltFit2Data_' xvar ap];
    assignin('base', plotName, figure)
    set(eval(plotName), 'Units', 'inches')
    set(eval(plotName), 'Position', [0 0 12 4])
    % ungrouped data
    subplot(nrows, ncols, 1)
    plot_fitToNutrient_sorted(xvar, Data, modData, ...
        'colDat', colDat, 'colMod', colMod, ...
        'standardised', standardised, 'connectors', connectors);
    % grouped by depth -- boxplot
    subplot(nrows, ncols, 2)
    plot_fitToNutrient_depth(xvar, Data, modData, ...
        'colDat', colDat, 'colMod', colMod, ...
        'standardised', standardised);
    % grouped by sample event -- boxplot
    subplot(nrows, ncols, 3)
    plot_fitToNutrient_event(xvar, Data, modData, ...
        'colDat', colDat, 'colMod', colMod, ...
        'standardised', standardised);
    switch standardised
        case false, Title = ['Model fit to ' xLab ' data'];
        case true, Title = ['Model fit to standardised ' xLab ' data'];
    end    
    sgtitle(Title)
end


switch save, case true
    % Group plots -- all scalar data
    plotName = 'pltFit2Data_raw';
    filename = 'fitToData_raw.png';
    if (exist(plotName, 'var') && isvalid(eval(plotName)))
        print(eval(plotName), fullfile(folder, filename), '-r300', '-dpng');
    end
    plotName = 'pltFit2Data_standardised';
    filename = 'fitToData_standardised.png';
    if (exist(plotName, 'var') && isvalid(eval(plotName)))
        print(eval(plotName), fullfile(folder, filename), '-r300', '-dpng');
    end
    % Individual plots -- separate scalar data types
    for i = 1:length(Vars)
        plotName = ['pltFit2Data_' Vars{i}];
        filename = ['fitToData_' Vars{i} '.png'];
        if (exist(plotName, 'var') && isvalid(eval(plotName)))
            print(eval(plotName), fullfile(folder, filename), '-r300', '-dpng');
        end
        plotName = ['pltFit2Data_' Vars{i} '_standardised'];
        filename = ['fitToData_' Vars{i} '_standardised.png'];
        if (exist(plotName, 'var') && isvalid(eval(plotName)))
            print(eval(plotName), fullfile(folder, filename), '-r300', '-dpng');
        end
    end
end


% Size data.

% Still need to include data variability in these plots...

% Fit to total biovolume within size class intervals (the data used to fit 
% the model)
trophicGroups = {'autotroph', 'heterotroph'}; pn = {'P','Z'};
xvar = 'BioVol';
logScale = 'semilogx'; % for size spectra data choose logScale = 'loglog' or 'semilogx'
waterOrigin = 'Atlantic';
% waterOrigin = 'Arctic';
connectors = true; % lines linking data to modelled equivalents
meanOnly = true; % display mean (over sample events) of modelled values -- for cleaner plot

for i = 1:length(trophicGroups)
    trophicGroup = trophicGroups{i};
    
    plotName = ['pltBioVol_Atl_' pn{i}];
    assignin('base', plotName, figure)
    set(eval(plotName), 'Units', 'inches')
    set(eval(plotName), 'Position', [0 0 16 6])

    subplot(1,3,1)
    plot_fitToSizeIntegrated_relative(xvar, trophicGroup, Data, modData, ...
        'logScale', logScale, 'waterOrigin', waterOrigin, ...
        'colDat', colDat, 'colMod', colMod,...
        'connectors', connectors, 'meanOnly', meanOnly, ...
        'includeLegend', true, 'legendLocation', 'northwest');
    
    subplot(1,3,2)
    plot_fitToSizeIntegrated_total(xvar, trophicGroup, Data, modData, ...
        'waterOrigin', waterOrigin, ...
        'colDat', colDat, 'colMod', colMod, ...
        'connectors', connectors, 'meanOnly', meanOnly);
    
    subplot(1,3,3)
    plot_fitToSizeIntegrated_absolute(xvar, trophicGroup, Data, modData, ...
        'logScale', logScale, 'waterOrigin', waterOrigin, ...
        'colDat', colDat, 'colMod', colMod, ...
        'connectors', connectors, 'meanOnly', meanOnly);
    
    switch xvar
        case 'BioVol', Xvar = 'biovolume';
        case 'CellConc', Xvar = 'cell concentration';
    end
    sgtitle(['Model fit to ' trophicGroup ' ' Xvar ' data'])
    
end


switch save, case true
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


close all

%~~~~~~~~~~~~~~~~~~~~
% Older plotting code...

% Summary plots displaying model fit to data
logPlot = true; % for scalar data choose logPlot = true or false
pltChl = plot_fitToData('chl_a', Data, modData, logPlot); pause(0.25)
pltPON = plot_fitToData('PON', Data, modData, logPlot); pause(0.25)
pltPOC = plot_fitToData('POC', Data, modData, logPlot); pause(0.25)
logPlot = false;
pltN = plot_fitToData('N', Data, modData, logPlot); pause(0.25)

% logPlot = 'loglog'; % for size spectra data choose logPlot = 'loglog' or 'semilogx'
logPlot = 'semilogx'; % for size spectra data choose logPlot = 'loglog' or 'semilogx'

% Comment out Arctic plots because we're fitting to Atlantic data
pltBioVol_Atl_P = plot_fitToData('BioVol', Data, modData, logPlot, ... 
    'trophicGroup', 'autotroph', 'waterOrigin', 'Atlantic'); pause(0.25)
% pltBioVol_Arc_P = plot_fitToData('BioVol', Data, modData, logPlot, ... 
%     'trophicGroup', 'autotroph', 'waterOrigin', 'Arctic'); pause(0.25)
pltBioVol_Atl_Z = plot_fitToData('BioVol', Data, modData, logPlot, ... 
    'trophicGroup', 'heterotroph', 'waterOrigin', 'Atlantic'); pause(0.25)
% pltBioVol_Arc_Z = plot_fitToData('BioVol', Data, modData, logPlot, ... 
%     'trophicGroup', 'heterotroph', 'waterOrigin', 'Arctic'); pause(0.25)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
    if exist('pltBioVol_Arc_P', 'var') && isvalid(pltBioVol_Arc_P)
        filename = 'fitToData_linePlots_BioVolSpectra_Arc_P.png';
        print(pltBioVol_Arc_P, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol_Arc_Z', 'var') && isvalid(pltBioVol_Arc_Z)
        filename = 'fitToData_linePlots_BioVolSpectra_Arc_Z.png';
        print(pltBioVol_Arc_Z, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol_Atl_P', 'var') && isvalid(pltBioVol_Atl_P)
        filename = 'fitToData_linePlots_BioVolSpectra_Atl_P.png';
        print(pltBioVol_Atl_P, fullfile(folder, filename), '-r300', '-dpng');
    end
    if exist('pltBioVol_Atl_Z', 'var') && isvalid(pltBioVol_Atl_Z)
        filename = 'fitToData_linePlots_BioVolSpectra_Atl_Z.png';
        print(pltBioVol_Atl_Z, fullfile(folder, filename), '-r300', '-dpng');
    end
    
end




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Residual errors -- put this code into a stand-alone function...

plt_scalarResid = figure;
plt_scalarResid.Units = 'inches';
plt_scalarResid.Position = [0 0 12 8];

colDat = [0.5, 0.5, 0.5];
deepLimit = 100;
shallowLimit = 20;
colDeep = [0, 0, 1];
colShallow = [0, 1, 0];
colHighlight = [1, 0, 0];
highlightDeep = Data.scalar.Depth >= deepLimit;
highlightShallow = Data.scalar.Depth<= shallowLimit;
% highlight = Data.scalar.Depth >= 100;
plotStandardisedData = true; % plot raw or standardised data
nv = length(Vars);

for ii = 1:nv
    subplot(2,2,ii)
    varLabel = Vars{ii};
    ind = strcmp(Data.scalar.Variable, varLabel);
    switch plotStandardisedData
        case true
            yobs = Data.scalar.scaled_Value(ind);
            ymod = modData.scalar.scaled_Value(ind,:);
            mainTitle = 'Residual errors: model minus standardised data';
        case false
            yobs = Data.scalar.Value(ind);
            ymod = modData.scalar.Value(ind,:);
            mainTitle = 'Residual errors: model minus data';
    end
    n = length(yobs);
    p = (1:n) ./ n;
    [yobs_sort, o] = sort(yobs);
    ymod_sort = ymod(o,:);
%     [ymod_sort, o] = sort(ymod);
%     yobs_sort = yobs(o,:);
    res = ymod_sort - yobs_sort;
    higlightAll = false(length(res),1);
    if (exist('highlightDeep', 'var') && any(highlightDeep(ind))) || ...
            (exist('highlightShallow', 'var') && any(highlightShallow(ind)))
        if (exist('highlightDeep', 'var') && any(highlightDeep(ind)))
            highlight_ = highlightDeep(ind);
            highlight_ = highlight_(o);
            higlightAll = higlightAll | highlight_;
            sDeep = scatter(yobs_sort(highlight_), res(highlight_), ...
                'MarkerEdgeColor', colDeep, 'MarkerFaceColor', colDeep);
            hold on
        end
        if exist('highlightShallow', 'var') && any(highlightShallow(ind))
            highlight_ = highlightShallow(ind);
            highlight_ = highlight_(o);
            higlightAll = higlightAll | highlight_;
            sShallow = scatter(yobs_sort(highlight_), res(highlight_), ...
                'MarkerEdgeColor', colShallow, 'MarkerFaceColor', colShallow);
            if ~(exist('highlightDeep', 'var') && any(highlightDeep(ind)))
                hold on
            end
        end
        if ii == 1
            leg = legend;
        end
    end
    scatter(yobs_sort(~higlightAll), res(~higlightAll), ... 
        'MarkerEdgeColor', colDat, 'LineWidth', 2)
%     scatter(1:n, res, 'MarkerEdgeColor', colDat, 'LineWidth', 2)
%     scatter(yobs_sort, res, 'MarkerEdgeColor', colDat, 'LineWidth', 2)
%     if exist('highlight', 'var') && any(highlight(ind))
%         highlight_ = highlight(ind);
%         highlight_ = highlight_(o);
%         scatter(yobs_sort(highlight_), res(highlight_), ... 
%             'MarkerEdgeColor', colHighlight, 'MarkerFaceColor', colHighlight)
%     end
    
    gc = gca;
    xl = gc.XLim;
    plot(xl, [0, 0], '--r')
    hold off
    xlabel(varLabel)
    ylabel('residual error')
    if ii == nv
        sgtitle(mainTitle)
    end
    if ii == 1
        % legend & text on panel
        if exist('sDeep', 'var') && exist('sShallow', 'var')
            legString = {['> ' num2str(deepLimit) ' m'], ... 
                ['< ' num2str(shallowLimit) ' m']};
            leg.String = legString;
        end
        if exist('sDeep', 'var') && ~exist('sShallow', 'var')
            legString = {['> ' num2str(deepLimit) ' m']};
            leg.String = legString;
        end
        if ~exist('sDeep', 'var') && exist('sShallow', 'var')
            legString = {['< ' num2str(shallowLimit) ' m']};
            leg.String = legString;
        end
        % include text
        yl = gc.YLim;
        xd = diff(gc.XLim);
        yd = diff(gc.YLim);
        text(xl(1) + 0.01 * xd, yl(2) - 0.05 * yd, 'overestimated', 'HorizontalAlignment', 'left')
        text(xl(1) + 0.01 * xd, yl(1) + 0.05 * yd, 'underestimated', 'HorizontalAlignment', 'left')
    end
end


% Size data - histograms .... see plots2.m
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


%% Map plots

% Particle trajectories from physical model

% Also see plots2.m

% ENSURE THAT THE FULL UNFILTERED FORCING DATA IS LOADED...

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




%% Contour plots -- single trajectories, or grouped by sample event

save = false;

% Choose one or more trajectories -- if multiple are selected then the plot
% will average over them.
sampleEvent = 1;
% All trajectories used for sampleEvent
traj = find(Data0.scalar.EventTraj(sampleEvent,:));
% If waterMass is either Atlantic OR Arctic then it may make sense to plot
% average over all trajectories, although there could be unwanted smoothing
% effects...
% Otherwise, if waterMass is a mixture of origins, then group trajectories
% by origin and make separate plots
waterMass = Data0.scalar.waterMass{sampleEvent};


switch waterMass
    case {'Atlantic', 'Arctic'}
        plt_Forc = plot_contour_DepthTime('forcing', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_DIN = plot_contour_DepthTime('DIN', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_OM = plot_contour_DepthTime('DOM_POM', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_N = plot_contour_DepthTime('phytoplankton_N', ...
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_Chl = plot_contour_DepthTime('phytoplankton_Chl', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_C = plot_contour_DepthTime('phytoplankton_C', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_N_C = plot_contour_DepthTime('phytoplankton_N_C', ...
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_P_Chl_N = plot_contour_DepthTime('phytoplankton_Chl_N', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_Z_C = plot_contour_DepthTime('zooplankton_C', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
            'Event', sampleEvent, 'waterOrigin', waterMass);
        plt_Z_N = plot_contour_DepthTime('zooplankton_N', ... 
            traj, out0, auxVars0, FixedParams, Forc0, 'linear', ...
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


close all



% Plot forcing data averaged over trajectories originating from Arctic & Atlantic

sampleEvent = 1;
% All trajectories used for sampleEvent
traj = find(Data0.scalar.EventTraj(sampleEvent,:));
% If waterMass is either Atlantic OR Arctic then it may make sense to plot
% average over all trajectories, although there could be unwanted smoothing
% effects...
% Otherwise, if waterMass is a mixture of origins, then group trajectories
% by origin and make separate plots
waterMass = Data0.scalar.waterMass{sampleEvent};


trajArctic = find(strcmp(Forc0.waterMass, 'Arctic'));
trajAtlantic = find(strcmp(Forc0.waterMass, 'Atlantic'));


plt_Forc = plot_contour_DepthTime('forcing', ...
    trajArctic, out0, auxVars0, FixedParams, Forc0, 'linear', ...
    'Event', sampleEvent, 'waterOrigin', 'Arctic');

plt_Forc = plot_contour_DepthTime('forcing', ...
    trajAtlantic, out0, auxVars0, FixedParams, Forc0, 'linear', ...
    'Event', sampleEvent, 'waterOrigin', 'Atlantic');





%% Time evolution

% These polygon plots should be extended to show water masses of Atlantic
% and of Arctic orgin because some sampling events use trajectories
% originating from both oceans... red for Atlantic blue for Arctic

% Choose event
sampleEvent = 1;
if ~ismember(sampleEvent, 1:Data0.scalar.nEvents), warning(['Choose event number within range (1, ' num2str(Data0.scalar.nEvents) ')']); end
% trajectory indices
traj = find(Data0.scalar.EventTraj(sampleEvent,:));

highlightColour = [1 0 1];
plotOptions = {'forcing', 'DIN', 'organicN', 'organicC', 'phytoplankton_C', ...
    'phytoplanktonStacked', 'phytoZooPlanktonStacked'};


% SOMETHING IS WRONG WITH THESE PLOTS... CODE NEEDS UPDATED...

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
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            % PAR
            subplot(3,1,3)
            plot_timeSeries_trajectoryPolygon('PAR', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'DIN'
            % DIN
            plt_DIN = figure;            
            plt_DIN.Units = 'inches';
            plt_DIN.Position = [0 0 12 8];
            
            subplot(2,1,1)
            plot_timeSeries_trajectoryPolygon('DIN', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,1,2)
            plot_timeSeries_trajectoryPolygon('DIN', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'organicN'
            plt_ON = figure;
            plt_ON.Units = 'inches';
            plt_ON.Position = [0 0 24 8];
            
            subplot(2,2,1)
            plot_timeSeries_trajectoryPolygon('DON', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,2)
            plot_timeSeries_trajectoryPolygon('PON', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,3)
            plot_timeSeries_trajectoryPolygon('DON', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,4)
            plot_timeSeries_trajectoryPolygon('PON', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'organicC'
            plt_OC = figure;
            plt_OC.Units = 'inches';
            plt_OC.Position = [0 0 24 8];

            subplot(2,2,1)
            plot_timeSeries_trajectoryPolygon('DOC', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,2)
            plot_timeSeries_trajectoryPolygon('POC', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'surface', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,3)
            plot_timeSeries_trajectoryPolygon('DOC', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            subplot(2,2,4)
            plot_timeSeries_trajectoryPolygon('POC', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'mean', 'highlightColour', highlightColour, 'plotNew', false);
            
        case 'phytoplankton_C'
            plt_P_C = plot_timeSeries_trajectoryPolygon('phytoplankton_C', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'max', 'highlightColour', highlightColour, 'fixedYaxis', false);
            plt_P_C_fixedScale = plot_timeSeries_trajectoryPolygon('phytoplankton_C', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0, ...
                'depth', 'max', 'highlightColour', highlightColour, 'fixedYaxis', true);
            
        case 'phytoplanktonStacked'
            plt_P_C_stacked = plot_timeSeries_trajectoryPolygon('phytoplanktonStacked', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0);

        case 'phytoZooPlanktonStacked'
            plt_P_C_stacked = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
                sampleEvent, traj, out0, auxVars0, FixedParams, Forc0, Data0);
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
% plt_Atlantic = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
%     [], [], out0, auxVars0, FixedParams, Forc0, Data0, ...
%     'waterMass', 'Atlantic', 'axesTextSize', axesTextSize, ... 
%     'titleTextSize', titleTextSize, 'legendTextSize', legendTextSize, ...
%     'legendPointSize', legendPointSize);
% plt_Atlantic.Units = 'inches';
% plt_Atlantic.Position = [0 0 10 4];
% gc = gca;
% yl = gc.YLim;
% 
% plt_Arctic = plot_timeSeries_trajectoryPolygon('phytoZooPlanktonStacked', ...
%     [], [], out0, auxVars0, FixedParams, Forc0, Data0, ...
%     'waterMass', 'Arctic', 'axesTextSize', axesTextSize, ... 
%     'titleTextSize', titleTextSize, 'legendTextSize', legendTextSize, ...
%     'legendPointSize', legendPointSize);
% plt_Arctic.Units = 'inches';
% plt_Arctic.Position = [0 0 10 4];
% gc = gca;
% gc.YLim = yl;

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


close all


%% Network plots -- fluxes, production

% Some of these network plots are too busy -- too many overlapping
% connections. The same information can be displayed differently as
% 'heatmaps' or tables...

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


% plt_OM = figure;
% plt_OM.Units = 'inches';
% plt_OM.Position = [0 0 10 3.8];
% subplot(1,2,1)
% plot_network('OMfluxes', 'carbon', auxVars0, FixedParams, Forc0, 'Arctic');
% subplot(1,2,2)
% plot_network('OMfluxes', 'carbon', auxVars0, FixedParams, Forc0, 'Atlantic');


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




% Try representing the feeding fluxes as a heatmap
Type = 'feedingFluxes';
nutrient = 'nitrogen';
traj = 'Atlantic';

plt_feedFlux_N = figure;
plt_feedFlux_N.Units = 'inches';
plt_feedFlux_N.Position = [0 0 10 7.5];

subplot(2,1,1)
plot_network_heatmap('feedingFluxes', 'nitrogen', auxVars0, FixedParams, Forc0, 'Atlantic');



subplot(2,1,2)
plot_network('feedingFluxes', 'nitrogen', auxVars0, FixedParams, Forc0, 'Atlantic');





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

