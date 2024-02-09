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
tag = '_RMS_Hellinger_ZPratio_Atlantic_final';
% tag = '_RMS_Hellinger2_Atlantic_aG_sigG_upweightAbnTot';
file_params = fullfile([fileName_params tag fileType_params]);

if ~exist(file_params, 'file'), warning('file_params does not exist -- check file names in results directory.'), file_params = []; end

% Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
%     'parFile', file_params);
Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', []);
display(Directories)

% Load saved outputs from optimisation runs or choose default set-up
loadFittedParams = true; % use output saved from optimisation run?
% Saved parameters file name and identifying tag
fileName_results = 'fittedParameters';
fileType_results = '.mat';
tag = '_RMS_Hellinger_ZPratio_Atlantic_final';
% tag = '_RMS_Hellinger2_Atlantic_aG_sigG_upweightAbnTot';
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

Params.Tref
Params.A
Params.h
Params.m2 = 0.05; % are code equations as they describe?
Params.aP
Params.theta
Params.xi
Params.aG
Params.k_G = 2.75;
Params.delta_opt
Params.sigG
Params.Lambda
Params.lambda_max
FixedParams.m_min
Params.rDOC = 0.04;
Params.rDON = 0.04;
Params.rPOC
Params.rPON
Params.beta1
Params.beta2
Params.beta3
Params.wPOM1 = 10;
Params.Qmin_QC_a = 0.035;
Params.Qmin_QC_b
Params.m_a
Params.m_b
Params.Qmax_delQ_a
Params.Qmax_delQ_b
Params.Vmax_QC_a
Params.Vmax_QC_b
Params.aN_QC_a
Params.aN_QC_b
Params.Gmax_a = 11;
Params.Gmax_b

% Params.Vmax_QC_a = 1e-9 / 14 * 0.024 / Params.Q_C_a
% Params.Vmax_QC_b = 1.1 - Params.Q_C_b

% Params.Qmin_QC_a = 1e-9 / 14 * 0.032 / Params.Q_C_a;


Params = updateParameters(Params, FixedParams, 'm2',0.05,'k_G',2.75,'rDOC',0.04,...
    'rDON',0.04,'wPOM1',10,'Qmin_QC_a',0.035,'Gmax_a',11, 'Vmax_QC_b', 1.1 - Params.Q_C_b, 'Qmin_QC_a', 0.1423);

% Initialise variables
if ~exist('v0', 'var') || ~isnumeric(v0)
    % If initial condition v0 has not been loaded then create initials.
    v0 = initialiseVariables(FixedParams, Params, Forc);
end
v00 = initialiseVariables(FixedParams, Params, Forc0); % create initials for full data (Forc0)

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Generate model outputs over trajectories linked to all data (Atlantic & Arctic) ...
clear out out0 auxVars auxVars0 modData0 modData
tic; disp('.. started at'); disp(datetime('now'))
[out0, auxVars0] = integrateTrajectories(FixedParams, Params, Forc0, v00, ... 
    FixedParams.odeIntegrator, FixedParams.odeOptions);
toc
% ... and only over the trajectories used to fit the model
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

savePlots = true;
folder = Directories.plotDir; % save plots here


%% Display data

% Nutrient & organic matter at depth
dataPlots.scalar = figure;
set(dataPlots.scalar, {'Units', 'Position'}, {'inches', [0 0 8 6]})

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


% Standardised data
dataPlots.scalarStandardised = figure;
set(dataPlots.scalarStandardised, {'Units', 'Position'}, {'inches', [0 0 6 9]})

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
plot_standardisedData2('scalar', 'N', Data0, 'covariate', 'Depth', ... 
    'pointAlpha', 0.4, 'densityCurve', true, ...
    'includeLegend', false, 'legendPosition', 'west');
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - ph, pw, ph])
plot_standardisedData2('scalar', 'N', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + (pwt - pw), 1 - pex - pht - ph, pw, ph])
plot_standardisedData2('scalar', 'PON', Data0, 'covariate', 'Depth', ... 
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - pht - ph, pw, ph])
plot_standardisedData2('scalar', 'PON', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + (pwt - pw), 1 - pex - 2 * pht - ph, pw, ph])
plot_standardisedData2('scalar', 'POC', Data0, 'covariate', 'Depth', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - 2 * pht - ph, pw, ph])
plot_standardisedData2('scalar', 'POC', Data0, 'covariate', 'Event', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + (pwt - pw), 1 - pex - 3 * pht - ph, pw, ph]);
plot_standardisedData2('scalar', 'chl_a', Data0, 'covariate', 'Depth', ...
    'pointAlpha', 0.4, 'densityCurve', true);
subplot('Position', [pex + 2 * pwt - pw, 1 - pex - 3 * pht - ph, pw, ph]);
plot_standardisedData2('scalar', 'chl_a', Data0, 'covariate', 'Event', ...
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


% Plot raw and standardised data in single figure
dataPlots.scalarAll = figure;
set(dataPlots.scalarAll, {'Units','Position'}, {'inches', [0 0 12 6]})

nrows = 2; % number of rows
ncols = 3; % and columns
pex = 0.02; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height (all panels and legend)
Pw = 1 - 2 * pex; % total plot width
yex = 0.2; % proportion of panel size devoted to axis labels -- y dimension
xex = 0.2;
pht = Ph / nrows; % panel height total
ph = pht - yex * pht; % panel height excluding labels
pwt = Pw / ncols;
pw = pwt - xex * pwt;

% Raw data vs depth
pointAlpha = 0.4;

subplot('Position', [pex + (pwt - pw), 1 - pex - ph, pw, ph])
plot_rawData('scalar', 'N', Data0, 'pointAlpha', pointAlpha, ...
    'includeLegend', true, 'legendTitle', 'water origin', 'legFontWeight', 'normal');

subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph, pw, ph])
plot_rawData('scalar', 'chl_a', Data0, 'pointAlpha', pointAlpha);

subplot('Position', [pex + (pwt - pw), 1 - pex - ph - pht, pw, ph])
plot_rawData('scalar', 'PON', Data0, 'pointAlpha', pointAlpha);

subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph - pht, pw, ph])
plot_rawData('scalar', 'POC', Data0, 'pointAlpha', pointAlpha);

% Standardised data vs event
group = {'N','PON','POC','chl_a'};
legLab = {'DIN','PON','POC','chl a'};
legPosition = 'northeast';
legOrientation = 'horizontal';
legColumns = 2;
legTitle = 'measurement';

subplot('Position', [pex + (3 * pwt - pw), 1 - pex - ph, pw, ph])
plot_standardisedData2('scalar', group, Data0, 'covariate', 'Depth', ... 
    'pointAlpha', pointAlpha, 'includeLegend', true, ... 
    'legendType', 'variables', 'legLab', legLab, 'legPosition', legPosition, ...
    'legOrientation', legOrientation, 'legColumns', legColumns, ...
    'legTitle', legTitle, 'legFontWeight', 'normal', 'logScale', true, ... 
    'YLim', [-350, -0.8]);

densityCurve = true;
lineAlpha = 0.5;
lineWidth = 1;
% legLab = {'standard normal', 'standardised data'};
legLab = {'N(0,1)', 'data'};
legTitle = 'distribution';
legPosition = 'northeast';
legOrientation = 'horizontal';
% legColumns = 2;
lineLim = [0, Data0.scalar.nEvents];

subplot('Position', [pex + (3 * pwt - pw), 1 - pex - ph - pht, pw, ph])
plot_standardisedData2('scalar', group, Data0, 'covariate', 'Event', ...
    'pointAlpha', pointAlpha, 'includeLegend', true, ...
    'legendType', 'distributions', 'legLab', legLab, 'legPosition', legPosition, ...
    'legTitle', legTitle, 'legOrientation', legOrientation, ... 
    'legFontWeight', 'normal', 'densityCurve', densityCurve, ... 
    'lineAlpha', lineAlpha, 'lineWidth', lineWidth, 'YLim', [0 48], 'lineLim', lineLim);


% Size spectra -- averaged over samples
dataPlots.sizeSpectra = figure;
set(dataPlots.sizeSpectra, {'Units', 'Position'}, {'inches', [0 0 8 6]})

% meanType = 'arithmetic';
meanType = 'geometric';
legendPosition = 'south';
subplot(2,1,1)
plot_rawData('sizeSpectra', 'CellConc', Data0, 'meanType', meanType);
subplot(2,1,2)
plot_rawData('sizeSpectra', 'BioVol', Data0, 'meanType', meanType, ...
    'includeLegend', true, 'legendPosition', legendPosition);


% Size spectra -- multipanel plots
dataPlots.sizeSpectraCellConc = figure;
set(dataPlots.sizeSpectraCellConc, {'Units', 'Position'}, {'inches', [0 0 16 12]})
Type = 'CellConc';

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


dataPlots.sizeSpectraBiovolume = figure;
set(dataPlots.sizeSpectraBiovolume, {'Units', 'Position'}, {'inches', [0 0 16 12]})
Type = 'BioVol';

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
dataPlots.sizeBinned = figure;
set(dataPlots.sizeBinned, {'Units', 'Position'}, {'inches', [0 0 8 6]})

subplot(2,1,1)
plot_rawData('sizeBinned', 'CellConc', Data0, 'pointAlpha', 0.5, ... 
    'includeLegend', true, 'legendPosition', 'south');
subplot(2,1,2)
plot_rawData('sizeBinned', 'BioVol', Data0, 'pointAlpha', 0.5);


switch savePlots, case true
    % save figures stored in dataPlots?
    if exist('dataPlots', 'var')
        fields = fieldnames(dataPlots);
        for i = 1:length(fields)
            if isvalid(dataPlots.(fields{i}))
                filename = ['data_' fields{i} '.png'];
                print(dataPlots.(fields{i}), fullfile(folder, filename), '-r300', '-dpng');
            end
        end
    end
end

close all
clear dataPlots


%% Model fit to data

%~~~~~~~~~~~~~~~~~~~~~~~
% Nutrient (scalar) data
%~~~~~~~~~~~~~~~~~~~~~~~

% Boxplots of variables vs either depth or sample event

fit2DataPlots.nutrient = figure;
set(fit2DataPlots.nutrient, {'Units', 'Position'}, {'inches', [0 0 16 8]})

standardised = true;
Vars = {'N','PON','POC','chl_a'};
colDat = [0, 0, 0];
colMod = [0, 1, 0];
% omitSingles = true; % exclude depths/events with only a single measurement?
omitSingles = 'merge'; % exclude depths/events with only a single measurement?

nrows = 2;
ncols = length(Vars);
% Position subplots to reduce excess white-space
pex = 0.01; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.2; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.15; % and y-axis
pht = Ph / nrows; % panel height total
ph = pht - yex * pht;
pwt = Pw / ncols;
pw = pwt - xex * pwt;

for j = 1:length(Vars)
    xvar = Vars{j};
    subplot('Position', [pex + (j * pwt - pw), 1 - pex - ph, pw, ph])
    plot_fitToNutrient_depth2(xvar, Data, modData, ...
        'colDat', colDat, 'colMod', colMod, ...
        'standardised', standardised, 'omitSingles', omitSingles)
    subplot('Position', [pex + (j * pwt - pw), 1 - pex - ph - pht, pw, ph])
    plot_fitToNutrient_event2(xvar, Data, modData, ...
        'colDat', colDat, 'colMod', colMod, ...
        'standardised', standardised, 'omitSingles', omitSingles)
end


%~~~~~~~~~~~~~~~~~~
% Size spectra data
%~~~~~~~~~~~~~~~~~~

fit2DataPlots.size = figure;
set(fit2DataPlots.size, {'Units', 'Position'}, {'inches', [0 0 13 8]})

nrows = 2;
ncols = 2;
% Position subplots to reduce excess white-space
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.125; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.125; % and y-axis
pht = Ph / nrows; % panel height total
ph = pht - yex * pht;
pwt = Pw / ncols;
pw = pwt - xex * pwt;

% depthLayer = mean([Data.size.dataBinned.DepthMin(1), Data.size.dataBinned.DepthMax(1)]);
depthLayer = Data.size.dataBinned.Depth(1);
% time = mean([Data.size.dataBinned.YeardayFirst(1), Data.size.dataBinned.YeardayLast(2)]);
time = Data.size.dataBinned.Yearday(1);
meanType = 'geometric';
compare2data = true;
normalised = true;
clear Cols
Cols.mod = [0 1 0];
Cols.obs = [0 0 0];
interp = true;
titleWeight = 'normal';
includeLegend = true;

subplot('Position', [pex + (pwt - pw), 1 - pex - ph, pw, ph])
xlab = [];
ylab = 'cell conc. density (log_{10}(ESD/1\mum)^{-1})';
Title = 'Autotrophs';
plotModelSizeSpectra2(out, auxVars, FixedParams, Data, Forc, 'spectrum',...
    'abundance', 'phyto', time, 'averaged', 'depthLayer', depthLayer, ...
    'meanType', meanType, 'compare2data', compare2data, ... 
    'normalised', normalised, 'ymin', 1e-10, 'xlab', xlab, 'ylab', ylab, ... 
    'Title', Title, 'titleWeight', titleWeight, 'Cols', Cols, ... 
    'interp', interp, 'includeLegend', includeLegend)

subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph, pw, ph])
ylab = '';
Title = 'Heterotrophs';
plotModelSizeSpectra2(out, auxVars, FixedParams, Data, Forc, 'spectrum',...
    'abundance', 'zoo', time, 'averaged', 'depthLayer', depthLayer, ...
    'meanType', meanType, 'compare2data', compare2data, ... 
    'normalised', normalised, 'ymin', 1e-10, 'xlab', xlab , 'ylab', ylab, ...
    'Title', Title, 'titleWeight', titleWeight, 'Cols', Cols, 'interp', interp)

subplot('Position', [pex + (pwt - pw), 1 - pex - pht - ph, pw, ph])
ylab = 'biovolume density (log_{10}(ESD/1\mum)^{-1})';
Title = [];
plotModelSizeSpectra2(out, auxVars, FixedParams, Data, Forc, 'spectrum',...
    'biovolume', 'phyto', time, 'averaged', 'depthLayer', depthLayer, ...
    'meanType', meanType, 'compare2data', compare2data, ... 
    'normalised', normalised, 'ymin', 1e-8, 'ylab', ylab, 'Title', Title, ...
    'Cols', Cols, 'interp', interp)

subplot('Position', [pex + (2 * pwt - pw), 1 - pex - pht - ph, pw, ph])
ylab = '';
plotModelSizeSpectra2(out, auxVars, FixedParams, Data, Forc, 'spectrum',...
    'biovolume', 'zoo', time, 'averaged', 'depthLayer', depthLayer, ...
    'meanType', meanType, 'compare2data', compare2data, ... 
    'normalised', normalised, 'ymin', 1e-8, 'ylab', ylab, 'Title', Title, ...
    'Cols', Cols, 'interp', interp)


switch savePlots, case true
    % save plots stored in fit2DataPlots
    if exist('fit2DataPlots', 'var')
        fields = fieldnames(fit2DataPlots);
        for i = 1:length(fields)
            if isvalid(fit2DataPlots.(fields{i}))
                filename = ['fit2Data_' fields{i} '.png'];
                print(fit2DataPlots.(fields{i}) , fullfile(folder, filename), '-r300', '-dpng');
            end
        end
    end
end


close all
clear fit2DataPlots



%% Fitted parameters

% Display fitted parameters in relation to their bounding values (in the
% table, columns widths shoukd be adjustable).

plt = plot_fittedParameters(results.optPar_summary);

switch savePlots, case true
    if exist('plt', 'var') && isvalid(plt)
        filename = 'fittedParameters.png';
        print(plt, fullfile(folder, filename), '-r300', '-dpng');
    end
end


close all
clear plt


%% Map plots

% generate struct containing all trajectories (unfiltered forcing data)
Forc_ = forcingSetUp(Directories, FixedParams, 'year', 2018);

plt_Map = figure;
set(plt_Map, {'Units', 'Position'}, {'inches', [0 0 8 8]})

% Plot all trajectories (Forc_) and highlight those selected for
% parameter-tuning (Forc0)
axes('position', [0.15, 0.15, 0.85, 0.85])
highlightTraj_ = nan(1,Forc0.nTraj);
x_ = Forc_.x;
y_ = Forc_.y;
for ii = 1:Forc0.nTraj
    x0 = Forc0.x(:,ii);
    y0 = Forc0.y(:,ii);    
    xi = find(all(x0 == x_));
    yi = find(all(y0 == y_));
    if xi == yi
        highlightTraj_(ii) = xi;
    end
end
highlightTraj = false(1,Forc_.nTraj);
highlightTraj(highlightTraj_) = true;

projection = 'lambert';
% alphaLine = 0.01;
alphaLine = 0.0075;
includeLegend = true;
legendPosition = 'east';
legendTextSize = 12;
legendTitle = 'Water origin';
legendTitleSize = 12;
polygonLineWidth = 3;
stripedBorder = false;
highlightStart = false;

plot_trajectoryMap(Directories, Forc_, 'projection', projection, ...
    'alphaLine', alphaLine, 'Data', Data0, 'newPlot', false, ...
    'includeLegend', includeLegend, 'legendPosition', legendPosition, ...
    'legendTitle', legendTitle, 'legendTitleSize', legendTitleSize, ...
    'legendTextSize', legendTextSize, 'polygonLineWidth', polygonLineWidth, ...
    'stripedBorder', stripedBorder, 'highlightStart', highlightStart, ... 
    'highlightTraj', highlightTraj);

% Plot the in-situ data
plt_Map_base = figure; % store base map
plt_Map_base = imshow(plt_Map);

figure
print(plt_Map_base)



axes('position', [0.1, 0.12, 0.5, 0.5])

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
    'latSpace', latSpace, 'Forc', Forc_, ... 
    'showWaterOrigin', showWaterOrigin, 'polygonAlpha', polygonAlpha, ...
    'polygonExpand', polygonExpand, 'polygonSmooth', polygonSmooth, ...
    'colourByDataType', colourByDataType, 'pieSize', pieSize, ...
    'includeLegend', includeLegend, 'legendPosition', legendPosition, ... 
    'omitMapGrid', omitMapGrid, 'colSat', colSat, ... 
    'fullAreaPolygon', fullAreaPolygon, 'polygonLineWidth', polygonLineWidth, ... 
    'trimPolygons', trimPolygons, 'legendTitle', legendTitle, ... 
    'legendTitleSize', legendTitleSize, 'legendTextSize', legendTextSize, ...
    'stripedBorder', stripedBorder, 'newPlot', false);

switch savePlots, case true
    if exist('plt_Map', 'var') && isvalid(plt_Map)
        filename = 'map_shipDataAndTrajectories.png';
        print(plt_Map, fullfile(folder, filename), '-r300', '-dpng');

        filename = 'map_shipDataAndTrajectories_res300.eps';
        exportgraphics(plt_Map, fullfile(folder, filename), 'Resolution', 300)

        filename = 'map_shipDataAndTrajectories.eps';
        exportgraphics(plt_Map, fullfile(folder, filename), 'Resolution', 3000)

    end
end


close all

%% Contour plots -- depth vs time

colourMap = 'plasma';
% colourMap = 'viridis';
% colourMap = 'magma';
% colourMap = 'inferno';
% colourMap = 'parula';

%~~~~~~~~~~~~~
% Forcing data 
%~~~~~~~~~~~~~

% Average over trajectories originating from Arctic & Atlantic.
% Use function optimisationOptions.m to filter data by region/water mass
fitTrajectories = 'Atlantic';
[~, ~, Forc_Atlantic, Data_Atlantic] = ...
    optimisationOptions(FixedParams, Params, Forc0, Data0, ...
    'fitTrajectories', fitTrajectories);
v0_Atlantic = initialiseVariables(FixedParams, Params, Forc_Atlantic);
Forc_Atlantic.integrateFullTrajectory = true;
% Generate model outputs over trajectories linked to Atlantic or Arctic
tic; disp('.. started at'); disp(datetime('now'))
[out_Atlantic, auxVars_Atlantic] = integrateTrajectories(FixedParams, Params, ... 
    Forc_Atlantic, v0_Atlantic, FixedParams.odeIntegrator, FixedParams.odeOptions);
toc
% repeat for Arctic
fitTrajectories = 'Arctic';
[~, ~, Forc_Arctic, Data_Arctic] = ...
    optimisationOptions(FixedParams, Params, Forc0, Data0, ...
    'fitTrajectories', fitTrajectories);
v0_Arctic = initialiseVariables(FixedParams, Params, Forc_Arctic);
Forc_Arctic.integrateFullTrajectory = true;
% Generate model outputs over trajectories linked to Atlantic or Arctic
tic; disp('.. started at'); disp(datetime('now'))
[out_Arctic, auxVars_Arctic] = integrateTrajectories(FixedParams, Params, ... 
    Forc_Arctic, v0_Arctic, FixedParams.odeIntegrator, FixedParams.odeOptions);
toc

% Create plot
contourPlots.Forc = figure;
set(contourPlots.Forc, {'Units', 'Position'}, {'inches', [0 0 10 9]})

% Left column = Arctic, right column = Atlantic
% Top row = temperature, middle row = PAR, bottom row = diffusivity

nrows = 3;
ncols = 2;
% Position subplots to reduce excess white-space
legw = (1 / 5) * (1 / ncols); % extra width for legends
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.125; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.125; % and y-axis
pht = Ph / nrows; % panel height total
ph = pht - yex * pht;
pwt = (Pw - legw) / ncols;
pw = pwt - xex * pwt;

ax(1) = subplot('Position', [pex + (pwt - pw), 1 - pex - ph, pw, ph]);
Title = 'Arctic';
ColourBar = true;
unitTemperature = 'temperature (\circC)';
xLabel = '';
plot_forcing_contour_DepthTime('T', Forc_Arctic, auxVars_Arctic, FixedParams, ... 
    'Title', Title, 'ColourBar', ColourBar, 'ColourBarLabel', unitTemperature, ...
    'xLabel', xLabel);
c1 = caxis;

ax(2) = subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph, pw, ph]);
Title = 'Atlantic';
ColourBar = true;
xLabel = '';
yLabel = '';
plot_forcing_contour_DepthTime('T', Forc_Atlantic, auxVars_Atlantic, FixedParams, ... 
    'Title', Title, 'ColourBar', ColourBar, 'ColourBarLabel', unitTemperature, ...
    'xLabel', xLabel, 'yLabel', yLabel);
c2 = caxis;

ax(3) = subplot('Position', [pex + (pwt - pw), 1 - pex - ph - pht, pw, ph]);
maxDepth = 50;
nYTicks = 6;
ColourBar = true;
scaleFactor = 1e-6; % rescale from muEin/d/m2 to Ein/d/m2
unitIrradiance = 'PAR (Ein day^{-1} m^{-2})';
xLabel = '';
plot_forcing_contour_DepthTime('PAR', Forc_Arctic, auxVars_Arctic, FixedParams, ... 
    'ColourBar', ColourBar, 'ColourBarLabel', unitIrradiance, ...
    'xLabel', xLabel, 'maxDepth', maxDepth, 'nYTicks', nYTicks, 'scaleFactor', scaleFactor);
c3 = caxis;

ax(4) = subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph - pht, pw, ph]);
nYTicks = 6;
ColourBar = true;
xLabel = '';
yLabel = '';
plot_forcing_contour_DepthTime('PAR', Forc_Atlantic, auxVars_Atlantic, FixedParams, ... 
    'ColourBar', ColourBar, 'ColourBarLabel', unitIrradiance, ...
    'xLabel', xLabel, 'yLabel', yLabel, 'maxDepth', maxDepth, 'nYTicks', nYTicks, 'scaleFactor', scaleFactor);
c4 = caxis;

ax(5) = subplot('Position', [pex + (pwt - pw), 1 - pex - ph - 2 * pht, pw, ph]);
logScale = true;
ColourBar = true;
unitDiffusivity = 'diffusivity (m^2 day^{-1})';
plot_forcing_contour_DepthTime('K', Forc_Arctic, auxVars_Arctic, FixedParams, ... 
    'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
    'logScale', logScale);
c5 = caxis;

ax(6) = subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph - 2 * pht, pw, ph]);
logScale = true;
ColourBar = true;
yLabel = '';
plot_forcing_contour_DepthTime('K', Forc_Atlantic, auxVars_Atlantic, FixedParams, ... 
    'ColourBar', ColourBar, 'ColourBarLabel', unitDiffusivity, ...
    'yLabel', yLabel, 'logScale', logScale);
c6 = caxis;

% Align color scales across rows
cT = [min([c1(1) c2(1)]), max([c1(2) c2(2)])];
cI = [min([c3(1) c4(1)]), max([c3(2) c4(2)])];
cK = [min([c5(1) c6(1)]), max([c5(2) c6(2)])];

caxis(ax(1), cT)
colorbar(ax(1), 'off')
pause(0.1)
caxis(ax(2), cT)
ax(2).Position(3) = ax(1).Position(3);
pause(0.1)

caxis(ax(3), cI)
colorbar(ax(3), 'off')
pause(0.1)
caxis(ax(4), cI)
ax(4).Position(3) = ax(3).Position(3);

caxis(ax(5), cK)
colorbar(ax(5), 'off')
caxis(ax(6), cK)
colorbar(ax(6), 'off')
cb = colorbar(ax(6), 'EastOutside');
Ticks = cb.Ticks;
ti = ismember(Ticks, -10:1:10);
Ticks = Ticks(ti);
TickLabels = arrayfun(@(z) num2str(10 .^ z), Ticks(:), 'UniformOutput', false);
set(cb, {'Ticks', 'TickLabels'}, {Ticks, TickLabels})
cb.Label.String = unitDiffusivity;
ax(6).Position(3) = ax(5).Position(3);

colormap(colourMap)


%~~~~~~~~~~~~~~
% Model outputs
%~~~~~~~~~~~~~~

% Display outputs for either Atlantic or Arctic waters -- averaged over all
% trajectories used for the selected water mass.

waterMass = 'Atlantic';

% DIN
contourPlots.DIN = figure;
set(contourPlots.DIN, {'Units', 'Position'}, {'inches', [0 0 8 3]})
var = 'DIN';
Title = [];
ColourBarLabel = 'DIN (mmol N m^{-3})';

plot_output_contour_DepthTime(var, out0 , auxVars0, FixedParams, Forc0, ...
    'waterOrigin', waterMass, 'Title', Title, 'ColourBarLabel', ColourBarLabel);
colormap(colourMap)


% DOM & POM
contourPlots.OM = figure;
set(contourPlots.OM, {'Units', 'Position'}, {'inches', [0 0 16 6]})
clear ax

nrows = 2;
ncols = 2;
% Position subplots to reduce excess white-space
legw = 0; % extra width for legends
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.125; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.125; % and y-axis
pht = Ph / nrows; % panel height total
ph = pht - yex * pht;
pwt = (Pw - legw) / ncols;
pw = pwt - xex * pwt;

ax(1) = subplot('Position', [pex + (pwt - pw), 1 - pex - ph, pw, ph]);
var = 'DOC';
Title = [];
ColourBarLabel = 'DOC (mmol C m^{-3})';
xlab = [];
plot_output_contour_DepthTime(var, out0 , auxVars0, FixedParams, Forc0, ...
    'waterOrigin', waterMass, 'Title', Title, 'ColourBarLabel', ColourBarLabel, ...
    'xlab', xlab);
c1 = caxis;

ax(2) = subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph, pw, ph]);
var = 'POC';
Title = [];
ColourBarLabel = 'POC (mmol C m^{-3})';
xlab = [];
ylab = [];
plot_output_contour_DepthTime(var, out0 , auxVars0, FixedParams, Forc0, ...
    'waterOrigin', waterMass, 'Title', Title, 'ColourBarLabel', ColourBarLabel, ...
    'xlab', xlab, 'ylab', ylab);
c2 = caxis;

ax(3) = subplot('Position', [pex + (pwt - pw), 1 - pex - ph - pht, pw, ph]);
var = 'DON';
Title = [];
ColourBarLabel = 'DON (mmol N m^{-3})';
plot_output_contour_DepthTime(var, out0 , auxVars0, FixedParams, Forc0, ...
    'waterOrigin', waterMass, 'Title', Title, 'ColourBarLabel', ColourBarLabel);
c3 = caxis;

ax(4) = subplot('Position', [pex + (2 * pwt - pw), 1 - pex - ph - pht, pw, ph]);
var = 'PON';
Title = [];
ColourBarLabel = 'PON (mmol N m^{-3})';
ylab = [];
plot_output_contour_DepthTime(var, out0 , auxVars0, FixedParams, Forc0, ...
    'waterOrigin', waterMass, 'Title', Title, 'ColourBarLabel', ColourBarLabel, ...
    'ylab', ylab);
c4 = caxis;

% Align color scales across rows
cC = [min([c1(1) c2(1)]), max([c1(2) c2(2)])];
cN = [min([c3(1) c4(1)]), max([c3(2) c4(2)])];

caxis(ax(1), cC)
pause(0.1)
caxis(ax(2), cC)
pause(0.1)

caxis(ax(3), cN)
pause(0.1)
caxis(ax(4), cN)

colormap(colourMap)
pause(0.1)


% Plankton -- separate plots for each trophic level and nutrient, and nutrient ratios
vars = {'P_N', 'P_Chl', 'P_C', 'P_biovolume', 'Z_N', 'Z_C', 'Z_biovolume', 'P_N_C', 'P_Chl_N', 'Z_N_C'};
smoothRatioColours = false;
xLim = [75, 300];
yLim = [1, 150];
Title = [];
xlab = 'year-day';
ylab = 'depth (m)';
sizeLab = true; % display size class within plot margins
sizeLabPosition = [0.02, 0.9, 0.2, 0.08]; % label position [x,y,width,height] specified as proportions of plot dimensions

% Subplot positions
legw = 0; % extra width for legends
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.125; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.125; % and y-axis

for i = 1:length(vars)
    var = vars{i};
    ratio = ismember(var, {'P_N_C', 'P_Chl_N', 'Z_N_C'});
    contourPlots.(var) = figure;    
    switch var(1)
        case 'P', nsize = FixedParams.nPP_size;
        case 'Z', nsize = FixedParams.nZP_size;
    end
    nc = floor(sqrt(nsize));
    nr = ceil(nsize / nc);
    set(contourPlots.(var), {'Units', 'Position'}, {'inches', [0, 0, nc * 8, nr * 3]})
    switch var
        case {'P_N', 'Z_N'}, ColourBarLabel = 'mmol N m^{-3}';
        case 'P_Chl', ColourBarLabel = 'mg Chl a m^{-3}';
        case {'P_C', 'Z_C'}, ColourBarLabel = 'mmol C m^{-3}';
        case {'P_N_C', 'Z_N_C'}, ColourBarLabel = 'mmol N (mmol C)^{-1}';
        case {'P_Chl_N'}, ColourBarLabel = 'mg Chl a (mmol N)^{-1}';
        case {'P_biovolume', 'Z_biovolume'}, ColourBarLabel = 'mm^3 m^{-3}';
    end
    
    % size label widths
    switch var(1)
        case 'P', ll = round([FixedParams.PPdia_intervals(1:end-1), FixedParams.PPdia_intervals(2:end)], 2, 'significant');
        case 'Z', ll = round([FixedParams.ZPdia_intervals(1:end-1), FixedParams.ZPdia_intervals(2:end)], 2, 'significant');
    end
    ll = arrayfun(@(z) num2str(z), ll, 'UniformOutput', false);
    ll = sum(strlength(ll), 2);
    ll = ll - min(ll);
    sizeLabWidths = sizeLabPosition(3) + 2*ll/100;

    clear ax cx
    
    pht = Ph / nr; % panel height total
    ph = pht - yex * pht;
    pwt = (Pw - legw) / nc;
    pw = pwt - xex * pwt;
    cx = cell(1,nsize);
    
    for sizeClass = 1:nsize
        jc = 1+mod(sizeClass-1,nc); % indexes column
        jj = repmat(1:nc, [nr, 1]);
        jr = jj(sizeClass); % indexes row
        ax(sizeClass) = subplot('Position', [pex + (jc * pwt - pw), 1 - pex - ph - (jr-1)*pht, pw, ph]);
        if jr == nr, xlab_ = xlab; else, xlab_ = []; end
        if jc == 1, ylab_ = ylab; else, ylab = []; end
        sizeLabPosition_ = sizeLabPosition;
        sizeLabPosition_(3) = sizeLabWidths(sizeClass);
        plot_output_contour_DepthTime(var, out0 , auxVars0, FixedParams, Forc0, ...
            'waterOrigin', waterMass, 'sizeClass', sizeClass, ...
            'ColourBarLabel', ColourBarLabel, 'xlab', xlab_, 'ylab', ylab_, ...
            'Title', Title, 'sizeLab', sizeLab, 'sizeLabPosition', sizeLabPosition_, ... 
            'xLim', xLim, 'yLim', yLim);
        cx{sizeClass} = caxis;
    end

    switch smoothRatioColours, case true
        % If var is a ratio then align color scales across subplots
        ci = [min(cellfun(@(z) z(1), cx)), max(cellfun(@(z) z(2), cx))];
        for j = 1:9
            caxis(ax(j), ci)
        end
    end
    
    colormap(colourMap)
    pause(0.1)
end


switch savePlots, case true
    % save plots stored in contourPlots
    if exist('contourPlots', 'var')
        fields = fieldnames(contourPlots);
        for i = 1:length(fields)
            if isvalid(contourPlots.(fields{i}))
                switch fields{i}
                    case 'Forc', filename = ['contourPlotDepthTime_' fields{i} '.png'];
                    otherwise, filename = ['contourPlotDepthTime_' fields{i} '_' waterMass '.png'];
                end
                print(contourPlots.(fields{i}), fullfile(folder, filename), '-r300', '-dpng');
            end
        end
    end
end


close all
clear contourPlots



%% Contour plots -- size vs time (abundance)

vars = {'P_N', 'P_Chl', 'P_C', 'Z_N', 'Z_C', 'P_cellDensity', 'Z_cellDensity', 'P_biovolume', 'Z_biovolume'};
waterMass = 'Atlantic';
colourMap = 'plasma';
logScale = false;
xLim = [100, 300];

for i = 1:length(vars)
    var = vars{i};
    trophicLevel = var(1);
    contourPlots.(var) = figure;
    set(contourPlots.(var), {'Units', 'Position'}, {'inches', [0, 0, 5, 3]})
    switch var
        case {'P_N', 'Z_N'}, ColourBarLabel = 'nitrogen (mmol N m^{-2})';
        case 'P_Chl', ColourBarLabel = 'chlorophyll a (mg chl a m^{-2})';
        case {'P_C', 'Z_C'}, ColourBarLabel = 'carbon (mmol C m^{-2})';
        case {'P_cellDensity', 'Z_cellDensity'}, ColourBarLabel = 'cell density (cells m^{-2})';
        case {'P_biovolume', 'Z_biovolume'}, ColourBarLabel = 'biovolume (mm^3 m^{-2})';
    end
    xlab = 'year-day';
    ylab = 'ESD (\mum)';
    switch trophicLevel
        case 'P', Title = 'autotrophs';
        case 'Z', Title = 'heterotrophs';
    end
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'logScale', logScale);
    colormap(eval(colourMap))
    pause(0.1)
end


% Group some plots and adjust colour scales
% Biovolumes
vars = {'P_biovolume', 'Z_biovolume'};
nvars = length(vars);
contourPlots.biovolumes = figure;
set(contourPlots.biovolumes, {'Units', 'Position'}, {'inches', [0, 0, 12, 4]})

logScale = true;
xLim = [75, 300];
% if ~logScale, minVal = 0; else, minVal = (1e-9*FixedParams.PPsize(1))*(1e3*FixedParams.Htot); end % 1 cell / L
minVal = 0;
minLevel = 0;

nlevels = 9; % contour levels
np = [1; 1]; % adjust contour levels -- default is np=1
lm = [4e4; 1.6e4]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
else
    levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
cp = 1 * ones(nvars, 1); % cp = 1 => default colours; other values rescale the colour map
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
colourMaps = cell(nvars,1);
for i = 1:nvars, colourMaps{i} = cm(Xnew(i,:),:); end

Title = [];

% Subplot positions
nr = 1;
nc = 2;
legw = 0.12;
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.1; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.1; % and y-axis
pht = Ph / nr; % panel height total
pwt = Pw / nc; % width
ph = pht - yex * pht;
pw = pwt - xex * pwt;
pwl = pw - legw * pwt;

clear ax
for i = 1:length(vars)
    ax(i) = subplot('Position', [pex + (i*pwt - pw), 1 - pex - ph, pwl, ph]);
    var = vars{i};
    trophicLevel = var(1);
    switch var
        case {'P_N', 'Z_N'}, ColourBarLabel = 'nitrogen (mmol N m^{-2})';
        case 'P_Chl', ColourBarLabel = 'chlorophyll a (mg chl a m^{-2})';
        case {'P_C', 'Z_C'}, ColourBarLabel = 'carbon (mmol C m^{-2})';
        case {'P_cellDensity', 'Z_cellDensity'}, ColourBarLabel = 'cell density (cells m^{-2})';
        case {'P_biovolume', 'Z_biovolume'}, ColourBarLabel = 'biovolume (mm^3 m^{-2})';
    end
    xlab = 'year-day';
    switch trophicLevel
        case 'P', ylab = 'autotroph ESD (\mum)';
        case 'Z', ylab = 'heterotroph ESD (\mum)';
    end
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'minVal', minVal, 'logScale', logScale, 'levels', levels(i,:));
    colormap(ax(i),colourMaps{i})
end


% Cell densities
vars = {'P_cellDensity', 'Z_cellDensity'};
nvars = length(vars);
contourPlots.cellDensities = figure;
set(contourPlots.cellDensities, {'Units', 'Position'}, {'inches', [0, 0, 12, 4]})

logScale = true;
xLim = [75, 300];
minVal = 0;
minLevel = 0;

nlevels = 9; % contour levels
np = [0.7; 0.7]; % adjust contour levels -- default is np=1
lm = [2.5e13; 2.5e11]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
else
    levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
cp = 1 * ones(nvars, 1); % cp = 1 => default colours; other values rescale the colour map
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
colourMaps = cell(nvars,1);
for i = 1:nvars, colourMaps{i} = cm(Xnew(i,:),:); end

Title = [];

clear ax
for i = 1:length(vars)
    ax(i) = subplot('Position', [pex + (i*pwt - pw), 1 - pex - ph, pwl, ph]);
    var = vars{i};
    trophicLevel = var(1);    
    switch var
        case {'P_N', 'Z_N'}, ColourBarLabel = 'nitrogen (mmol N m^{-2})';
        case 'P_Chl', ColourBarLabel = 'chlorophyll a (mg chl a m^{-2})';
        case {'P_C', 'Z_C'}, ColourBarLabel = 'carbon (mmol C m^{-2})';
        case {'P_cellDensity', 'Z_cellDensity'}, ColourBarLabel = 'cell density (cells m^{-2})';
        case {'P_biovolume', 'Z_biovolume'}, ColourBarLabel = 'biovolume (mm^3 m^{-2})';
    end    
    xlab = 'year-day';
    switch trophicLevel
        case 'P', ylab = 'autotroph ESD (\mum)';
        case 'Z', ylab = 'heterotroph ESD (\mum)';
    end
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ...
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'minVal', minVal, 'logScale', logScale, 'levels', levels(i,:));
    colormap(ax(i),colourMaps{i})
end


switch savePlots, case true
    if exist('contourPlots', 'var')
        fields = fieldnames(contourPlots);
        for i = 1:length(fields)            
            filename = ['contourPlotSizeTime_' fields{i} '_' waterMass '.png'];
            print(contourPlots.(fields{i}) , fullfile(folder, filename), '-r300', '-dpng');
        end
    end
end

close all
clear contourPlots


%% Contour plots -- size vs time (OM production)

vars = {'P_ON', 'P_OC', 'Z_ON', 'Z_OC'};
waterMass = 'Atlantic';
colourMap = 'plasma';
logScale = false;
xLim = [100, 300];
Title = [];

for i = 1:length(vars)
    var = vars{i};
    trophicLevel = var(1);
    contourPlots.(var) = figure;
    set(contourPlots.(var), {'Units', 'Position'}, {'inches', [0, 0, 5, 3]})
    switch var
        case {'P_ON', 'Z_ON'}, ColourBarLabel = 'organic nitrogen (mmol N m^{-2} d^{-1})';
        case {'P_OC', 'Z_OC'}, ColourBarLabel = 'organic carbon (mmol C m^{-2} d^{-1})';
    end
    xlab = 'year-day';
    switch trophicLevel
        case 'P', ylab = 'autotroph ESD (\mum)';
        case 'Z', ylab = 'heterotroph ESD (\mum)';
    end
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'logScale', logScale);
    colormap(eval(colourMap))
    pause(0.1)
end


% Group plots and adjust colour scales
% Total organic carbon & nitrogen
vars = {'P_OC', 'P_ON', 'Z_OC', 'Z_ON'};
nvars = length(vars);
contourPlots.OM = figure;
set(contourPlots.OM, {'Units', 'Position'}, {'inches', [0, 0, 12, 8]})

logScale = true;
xLim = [75, 300];
minVal = zeros(nvars, 1);
minLevel = 0;

nlevels = 9; % contour levels
np = [1; 1; 1; 1]; % adjust contour levels -- default is np=1
lm = [30; 7; 90; 14]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
else
    levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
cp = ones(nvars,1);
% cp = [1.5; 1.5; 1.5; 1.5]; % cp = 1 => default colours; other values rescale the colour map
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
for jj = 1:size(Xnew,1), colourMaps{jj} = cm(Xnew(jj,:),:); end

Title = [];

% Subplot positions
nr = 2;
nc = 2;
legw = 0.05;
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.1; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.1; % and y-axis
pht = Ph / nr; % panel height total
pwt = Pw / nc; % width
ph = pht - yex * pht;
pw = pwt - xex * pwt;
pwl = pw - legw * pwt;

clear ax
for i = 1:length(vars)
    jj = repmat(1:nc, [nr 1]);
    j = jj(i);
    ax(i) = subplot('Position', [pex + ((1 + mod(i-1,nc))*pwt - pw), 1 - pex - ph - (j-1)*pht, pwl, ph]);
    var = vars{i};
    trophicLevel = var(1);    
    switch var
        case {'P_ON', 'Z_ON'}, ColourBarLabel = 'organic nitrogen (mmol N m^{-2} d^{-1})';
        case {'P_OC', 'Z_OC'}, ColourBarLabel = 'organic carbon (mmol C m^{-2} d^{-1})';
    end
    if j == nr
        xlab = 'year-day';
    else
        xlab = [];
    end
    switch trophicLevel
        case 'P', ylab = 'autotroph ESD (\mum)';
        case 'Z', ylab = 'heterotroph ESD (\mum)';
    end
    if mod(i-1,nc) ~= 0
        ylab = [];
    end
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'minVal', minVal(i), 'logScale', logScale, 'levels', levels(i,:));
    
    colormap(ax(i),colourMaps{i})
end


% Separate DOM and POM, mortality and messy feeding
% Columns for DOM & POM, rows for trophic level/mortality & messy feeding 
vars = {'P_DOC', 'P_POC', 'Z_DOC_mort', 'Z_POC_mort', 'Z_DOC_mess', 'Z_POC_mess'};
nvars = length(vars);
contourPlots.OCdetailed = figure;
set(contourPlots.OCdetailed, {'Units', 'Position'}, {'inches', [0, 0, 12, 12]})

labels = [repmat({'mortality'}, [4,1]); repmat({'messy feeding'}, [2,1])];
labelCol = [1 1 1];

logScale = true;
xLim = [75, 300];
minVal = zeros(nvars, 1);
minLevel = 0;

nlevels = 9; % contour levels
np = [1; 1; 1; 1; 1; 1]; % adjust contour levels -- default is np=1
lm = [25; 5.5; 4.5; 3; 60; 22]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
else
    levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
% cp = [1.5; 1.5; 1.5; 1.5]; % cp = 1 => default colours; other values rescale the colour map
cp = ones(nvars, 1);
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
for jj = 1:size(Xnew,1), colourMaps{jj} = cm(Xnew(jj,:),:); end

Title = [];

% Subplot positions
nr = 3;
nc = 2;
legw = 0.05;
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.1; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.1; % and y-axis
pht = Ph / nr; % panel height total
pwt = Pw / nc; % width
ph = pht - yex * pht;
pw = pwt - xex * pwt;
pwl = pw - legw * pwt;

clear ax
for i = 1:nvars
    jc = 1 + mod(i-1, nc); % indexes column
    jj = repelem(1:nr, nc);
    jr = jj(i); % indexes row    
    ax(i) = subplot('Position', [pex + (jc*pwt - pw), 1 - pex - ph - (jr-1)*pht, pwl, ph]);
    var = vars{i};
    trophicLevel = var(1);        
    switch var
        case {'P_DON', 'Z_DON_mess', 'Z_DON_mort'}, ColourBarLabel = 'DON (mmol N m^{-2} d^{-1})';
        case {'P_PON', 'Z_PON_mess', 'Z_PON_mort'}, ColourBarLabel = 'PON (mmol N m^{-2} d^{-1})';
        case {'P_DOC', 'Z_DOC_mess', 'Z_DOC_mort'}, ColourBarLabel = 'DOC (mmol C m^{-2} d^{-1})';
        case {'P_POC', 'Z_POC_mess', 'Z_POC_mort'}, ColourBarLabel = 'POC (mmol C m^{-2} d^{-1})';
    end
    if jr == nr
        xlab = 'year-day';
    else
        xlab = [];
    end
    switch trophicLevel
        case 'P', ylab = 'autotroph ESD (\mum)';
        case 'Z', ylab = 'heterotroph ESD (\mum)';
    end
    if jc ~= 1
        ylab = [];
    end
    
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'minVal', minVal(i), 'logScale', logScale, 'levels', levels(i,:));
    
    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    yl = log10(yl);
    text(xl(1) + 0.05 * diff(xl), 10 .^ (yl(2) - 0.05 * diff(yl)), labels{i}, 'Color', labelCol)
    
    colormap(ax(i),colourMaps{i})

end


% Repeat for nitrogen...
vars = {'P_DON', 'P_PON', 'Z_DON_mort', 'Z_PON_mort', 'Z_DON_mess', 'Z_PON_mess'};
nvars = length(vars);
contourPlots.ONdetailed = figure;
set(contourPlots.ONdetailed, {'Units', 'Position'}, {'inches', [0, 0, 12, 12]})

labels = [repmat({'mortality'}, [4,1]); repmat({'messy feeding'}, [2,1])];

logScale = true;
xLim = [75, 300];
minVal = zeros(nvars, 1);
minLevel = 0;

nlevels = 9; % contour levels
np = [1; 1; 1; 1; 1; 1]; % adjust contour levels -- default is np=1
lm = [6; 2.25; 2; 1.6; 8; 3]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
else
    levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
% cp = [1.5; 1.5; 1.5; 1.5]; % cp = 1 => default colours; other values rescale the colour map
cp = ones(nvars, 1);
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
for jj = 1:size(Xnew,1), colourMaps{jj} = cm(Xnew(jj,:),:); end

Title = [];

% Subplot positions
nr = 3;
nc = 2;
legw = 0.05;
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.1; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.1; % and y-axis
pht = Ph / nr; % panel height total
pwt = Pw / nc; % width
ph = pht - yex * pht;
pw = pwt - xex * pwt;
pwl = pw - legw * pwt;

clear ax
for i = 1:nvars
    jc = 1 + mod(i-1, nc); % indexes column
    jj = repelem(1:nr, nc);
    jr = jj(i); % indexes row    
    ax(i) = subplot('Position', [pex + (jc*pwt - pw), 1 - pex - ph - (jr-1)*pht, pwl, ph]);
    var = vars{i};
    trophicLevel = var(1);        
    switch var
        case {'P_DON', 'Z_DON_mess', 'Z_DON_mort'}, ColourBarLabel = 'DON (mmol N m^{-2} d^{-1})';
        case {'P_PON', 'Z_PON_mess', 'Z_PON_mort'}, ColourBarLabel = 'PON (mmol N m^{-2} d^{-1})';
        case {'P_DOC', 'Z_DOC_mess', 'Z_DOC_mort'}, ColourBarLabel = 'DOC (mmol C m^{-2} d^{-1})';
        case {'P_POC', 'Z_POC_mess', 'Z_POC_mort'}, ColourBarLabel = 'POC (mmol C m^{-2} d^{-1})';
    end
    if jr == nr
        xlab = 'year-day';
    else
        xlab = [];
    end
    switch trophicLevel
        case 'P', ylab = 'autotroph ESD (\mum)';
        case 'Z', ylab = 'heterotroph ESD (\mum)';
    end
    if jc ~= 1
        ylab = [];
    end
%     plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
%         'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
%         'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
%         'minVal', minVal(i), 'logScale', logScale);
    plot_output_contour_SizeTime(var, out0 , auxVars0, FixedParams, Forc0, ...
        'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ... 
        'xlab', xlab, 'ylab', ylab, 'Title', Title, 'xLim', xLim, ...
        'minVal', minVal(i), 'logScale', logScale, 'levels', levels(i,:));
    
    xl = get(gca, 'XLim');
    yl = get(gca, 'YLim');
    yl = log10(yl);
    text(xl(1) + 0.05 * diff(xl), 10 .^ (yl(2) - 0.05 * diff(yl)), labels{i}, 'Color', labelCol)
    
    colormap(ax(i),colourMaps{i})

end


switch savePlots, case true
    if exist('contourPlots', 'var')
        fields = fieldnames(contourPlots);
        for i = 1:length(fields)            
            filename = ['contourPlotSizeTime_' fields{i} '_' waterMass '.png'];
            print(contourPlots.(fields{i}) , fullfile(folder, filename), '-r300', '-dpng');
        end
    end
end

close all
clear contourPlots


%% Contour plots -- size vs size (grazing)

% Can plot total grazing losses across full year (mmol x / m^2 / year), or
% select a time period (growing season) then plot daily average grazing
% losses within that period (mmol x / m^2 / day).

nutrients = FixedParams.ZP_nut;
vars = {'Z_P_N', 'Z_P_C', 'Z_Z_N', 'Z_Z_C'}; % pred_prey_element
waterMass = 'Atlantic';
colourMap = 'plasma';
% powerScale = 0.25; % plot rescaled values (x ^ powerScale) for clarity (smoothing the values with powerScale = 0.25 improves linear interpolation accuracy => better contours)
powerScale = 1; % plot rescaled values (x ^ powerScale) for clarity (smoothing the values with powerScale = 0.25 improves linear interpolation accuracy => better contours)
logScale = false;
Times = 100:250; % time period to plot
totOrAve = 'average'; % 'average' => average daily losses over Times, 'total' => total losses over Times
xlab = 'prey ESD (\mum)';
ylab = 'predator ESD (\mum)';
Title = [];

% Subplot positions
nr = 1;
nc = 2;
legw = 0.125;
pex = 0.025; % proportion of total plot size used for outer margins
Ph = 1 - 2 * pex; % total plot height
Pw = 1 - 2 * pex; % and width
xex = 0.125; % proportion of panel size deveoted to axis labels, x-axis
yex = 0.125; % and y-axis
pht = Ph / nr; % panel height total
pwt = Pw / nc; % width
ph = pht - yex * pht;
pw = pwt - xex * pwt;
pwl = pw - legw * pwt;

for j = 1:length(nutrients)
    nut = nutrients{j};
    switch nut
        case 'N', ColourBarLabel = 'mmol N m^{-2}';
        case 'C', ColourBarLabel = 'mmol C m^{-2}';
    end
    switch totOrAve
        case 'total', ColourBarLabel = [ColourBarLabel ' y^{-1}'];
        case 'average', ColourBarLabel = [ColourBarLabel ' d^{-1}'];
    end
    plotName = ['grazingLosses_' nut];
    contourPlots.(plotName) = figure;
    set(contourPlots.(plotName), {'Units', 'Position'}, {'inches', [0, 0, 10, 4]})
    vars_ = vars(contains(vars, nut));
    for i = 1:length(vars_)
        subplot('Position', [pex + (i*pwt - pw), 1 - pex - ph, pwl, ph]);
        var = vars_{i};
        prey = var(3);
        switch prey
            case 'P', xlab_ = ['autotroph ' xlab];
            case 'Z', xlab_ = ['heterotroph ' xlab];
        end
        if ~isempty(Title)
            switch prey
                case 'P', Title = 'autotroph grazing losses';
                case 'Z', Title = 'heterotroph grazing losses';
            end
        end
        plot_output_contour_SizeSize(var, out0 , auxVars0, FixedParams, Forc0, ...
            'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ...
            'xlab', xlab_, 'ylab', ylab, 'Title', Title, ...
            'logScale', logScale, 'Times', Times, 'totOrAve', totOrAve);
        
        colormap(colourMap)
    end
end


% Adjust scalings for grazing carbon losses and replot
vars = {'Z_P_C', 'Z_Z_C'}; % pred_prey_element
nvars = length(vars);
nutrients = {'C'};
waterMass = 'Atlantic';

logScale = true;
Times = 100:250; % time period to plot
totOrAve = 'average'; % 'average' => average daily losses over Times, 'total' => total losses over Times
xlab = 'prey ESD (\mum)';
ylab = 'predator ESD (\mum)';
Title = [];

nlevels = 9; % contour levels
minLevel = 0;
np = 1*ones(nvars, 1); % adjust contour levels -- default is np=1
lm = [11; 2.5]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
% else
%     levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
cp = 1*ones(nvars,1); % cp = 1 => default colours; other values rescale the colour map
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
colourMaps{1} = cm(Xnew(1,:),:);
colourMaps{2} = cm(Xnew(2,:),:);
colBarMax = [0.8, 0.5]; % controls max value displayed in colourbar

for j = 1:length(nutrients)
    nut = nutrients{j};
    switch nut
        case 'N', ColourBarLabel = 'mmol N m^{-2}';
        case 'C', ColourBarLabel = 'mmol C m^{-2}';
    end
    switch totOrAve
        case 'total', ColourBarLabel = [ColourBarLabel ' y^{-1}'];
        case 'average', ColourBarLabel = [ColourBarLabel ' d^{-1}'];
    end
    plotName = ['grazingLosses_' nut];
    contourPlots.(plotName) = figure;
    set(contourPlots.(plotName), {'Units', 'Position'}, {'inches', [0, 0, 10, 4]})
    vars_ = vars(contains(vars, nut));
    for i = 1:length(vars_)
        subplot('Position', [pex + (i*pwt - pw), 1 - pex - ph, pwl, ph]);
        var = vars_{i};
        prey = var(3);
        switch prey
            case 'P', xlab_ = ['autotroph ' xlab];
            case 'Z', xlab_ = ['heterotroph ' xlab];
        end
        if ~isempty(Title)
            switch prey
                case 'P', Title = 'autotroph grazing losses';
                case 'Z', Title = 'heterotroph grazing losses';
            end
        end
        plot_output_contour_SizeSize(var, out0 , auxVars0, FixedParams, Forc0, ...
            'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ...
            'xlab', xlab_, 'ylab', ylab, 'Title', Title, ...
            'logScale', logScale, 'Times', Times, 'totOrAve', totOrAve, ...
            'levels', levels(i,:), 'colBarMax', colBarMax(i));

        colormap(colourMaps{i})
    end
end


% Adjust scalings for grazing nitrogen losses and replot
vars = {'Z_P_N', 'Z_Z_N'}; % pred_prey_element
nvars=length(vars);
nutrients = {'N'};
waterMass = 'Atlantic';

logScale = true;
Times = 100:250; % time period to plot
totOrAve = 'average'; % 'average' => average daily losses over Times, 'total' => total losses over Times
xlab = 'prey ESD (\mum)';
ylab = 'predator ESD (\mum)';
Title = [];

nlevels = 9; % contour levels
minLevel = 0;
np = 1*ones(nvars, 1); % adjust contour levels -- default is np=1
lm = [2.5; 1.1]; % upper limits for contours
if logScale
    lm = log10(lm+1);
end
levels = minLevel + (lm - minLevel) .* linspace(0, 1, nlevels) .^ np;
if ~logScale
    levels = round(levels, 2, 'significant');
% else
%     levels = log10(round(10 .^ levels, 2, 'significant'));
end

% rescale colour map
cp = 1 * ones(nvars,1); % cp = 1 => default colours; other values rescale the colour map
cm = eval(colourMap);
cn = size(cm,1);
X = 1:cn;
Xnew = linspace(0, 1, length(X));
Xnew = round(cn .* Xnew .^ cp);
Xnew(Xnew == 0) = 1;
colourMaps{1} = cm(Xnew(1,:),:);
colourMaps{2} = cm(Xnew(2,:),:);
colBarMax = [0, 0]; % controls max value displayed in colourbar

for j = 1:length(nutrients)
    nut = nutrients{j};
    switch nut
        case 'N', ColourBarLabel = 'mmol N m^{-2}';
        case 'C', ColourBarLabel = 'mmol C m^{-2}';
    end
    switch totOrAve
        case 'total', ColourBarLabel = [ColourBarLabel ' y^{-1}'];
        case 'average', ColourBarLabel = [ColourBarLabel ' d^{-1}'];
    end
    plotName = ['grazingLosses_' nut];
    contourPlots.(plotName) = figure;
    set(contourPlots.(plotName), {'Units', 'Position'}, {'inches', [0, 0, 10, 4]})
    vars_ = vars(contains(vars, nut));
    for i = 1:length(vars_)
        subplot('Position', [pex + (i*pwt - pw), 1 - pex - ph, pwl, ph]);
        var = vars_{i};
        prey = var(3);
        switch prey
            case 'P', xlab_ = ['autotroph ' xlab];
            case 'Z', xlab_ = ['heterotroph ' xlab];
        end
        if ~isempty(Title)
            switch prey
                case 'P', Title = 'autotroph grazing losses';
                case 'Z', Title = 'heterotroph grazing losses';
            end
        end
        plot_output_contour_SizeSize(var, out0 , auxVars0, FixedParams, Forc0, ...
            'waterOrigin', waterMass, 'ColourBarLabel', ColourBarLabel, ...
            'xlab', xlab_, 'ylab', ylab, 'Title', Title, ...
            'logScale', logScale, 'Times', Times, 'totOrAve', totOrAve, ...
            'levels', levels(i,:), 'colBarMax', colBarMax(i));
        
        colormap(colourMaps{i})
    end
end


switch savePlots, case true
    if exist('contourPlots', 'var')
        fields = fieldnames(contourPlots);
        for i = 1:length(fields)            
            filename = ['contourPlotSizeSize_' fields{i} '_' waterMass '.png'];
            print(contourPlots.(fields{i}) , fullfile(folder, filename), '-r300', '-dpng');
        end
    end
end

close all
clear contourPlots



%% Time evolution

waterMass = 'Atlantic';

% Forcing data
timeSeriesPlots.forcing = figure;
set(timeSeriesPlots.forcing, {'Units', 'Position'}, {'inches', [0 0 8 12]})

polygonColour = [0.5 0.5 0.5];
lineColour = [0 0 0];

% temperature
subplot(3,1,1)
xlab = [];
plot_forcing_timeSeriesPolygon('temperature', FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour, 'xlab', xlab);
% diffusivity
subplot(3,1,2)
plot_forcing_timeSeriesPolygon('diffusivity', FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour, 'xlab', xlab);
% PAR
subplot(3,1,3)
plot_forcing_timeSeriesPolygon('PAR', FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'surface', 'polygonColour', polygonColour, ...
    'lineColour', lineColour);



% Inorganic nutrient
timeSeriesPlots.DIN = figure;
set(timeSeriesPlots.DIN, {'Units', 'Position'}, {'inches', [0 0 8 4]})

plot_output_timeSeriesPolygon('DIN', out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour);


% Organic nutrient
timeSeriesPlots.OM = figure;
set(timeSeriesPlots.OM, {'Units', 'Position'}, {'inches', [0 0 16 8]})

subplot(2,2,1)
var = 'DON';
plot_output_timeSeriesPolygon(var, out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour);

subplot(2,2,2)
var = 'DOC';
plot_output_timeSeriesPolygon(var, out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour);

subplot(2,2,3)
var = 'PON';
plot_output_timeSeriesPolygon(var, out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour);

subplot(2,2,4)
var = 'POC';
plot_output_timeSeriesPolygon(var, out0, auxVars0, FixedParams, Forc0, Data0, ...
    'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
    'lineColour', lineColour);


% Plankton -- separate plot for each trophic level & nutrient type

vars = {'P_N', 'P_Chl', 'P_C', 'Z_N', 'Z_C', 'P_N_C', 'P_Chl_N', 'Z_N_C'};
fixedYaxis = false; % same y-axis scale on all subplots?

sizeLab = true; % display size class within plot margins
sizeLabPosition = [0.02, 0.88, 0.22, 0.1]; % label position [x,y,width,height] specified as proportions of plot dimensions
sizeLabWidths = [0.22, 0.22, 0.25, 0.24, 0.23, 0.23, 0.23, 0.25, 0.27]; % manually set label widths... there must be a good method to automate this...

for i = 1:length(vars)
    var = vars{i};
    timeSeriesPlots.(var) = figure;    
    switch var(1)
        case 'P', nsize = FixedParams.nPP_size;
        case 'Z', nsize = FixedParams.nZP_size;
    end
    nc = floor(sqrt(nsize));
    nr = ceil(nsize / nc);
    set(timeSeriesPlots.(var), {'Units', 'Position'}, {'inches', [0, 0, nc * 8, nr * 4]})
    ratioPlot = ismember(var, {'P_N_C', 'P_Chl_N', 'Z_N_C'});
    YLim = nan(nsize,2);    
    for j = 1:nsize
        subplot(nr, nc, j)
        sizeLabPosition(3) = sizeLabWidths(j);
        plot_output_timeSeriesPolygon(var, out0, auxVars0, FixedParams, Forc0, Data0, ...
            'waterMass', waterMass, 'depth', 'mean', 'polygonColour', polygonColour, ...
            'lineColour', lineColour, 'sizeClass', j, 'sizeLab', sizeLab, ...
            'sizeLabPosition', sizeLabPosition);
        ax(j) = gca;
        YLim(j,:) = get(gca, 'YLim');
    end
    switch ratioPlot
        case false
            switch fixedYaxis
                case true, YLim = repmat([0 max(YLim(:,2))], [nsize 1]);
                case false, YLim = [zeros(nsize, 1) YLim(:,2)];
            end
    end
    for j = 1:nsize
        ax(j).YLim = YLim(j,:);
    end
    clear ax
end


% Similar plots that combine size classes by stacking polygons

waterMass = 'Atlantic';
var = 'PZ_C_stacked';

colsP = [0 100 0; 127 255 212] ./ 255; % dark green to aquamarine
colsZ = [128 0 0; 255 255 0] ./ 255;   % maroon to yellow

timeSeriesPlots.PZ_C_stacked = figure;
set(timeSeriesPlots.PZ_C_stacked, {'Units', 'Position'}, {'inches', [0 0 8 4]})

plot_output_timeSeriesPolygon(var, out0, auxVars0, ...
    FixedParams, Forc0, Data0, 'waterMass', waterMass, ...
    'colsP', colsP, 'colsZ', colsZ);


timeSeriesPlots.PZ_N_stacked = figure;
set(timeSeriesPlots.PZ_N_stacked, {'Units', 'Position'}, {'inches', [0 0 8 4]})
var = 'PZ_N_stacked';

plot_output_timeSeriesPolygon(var, out0, auxVars0, ...
    FixedParams, Forc0, Data0, 'waterMass', waterMass, ...
    'colsP', colsP, 'colsZ', colsZ);


timeSeriesPlots.PZ_biomass_stacked = figure;
set(timeSeriesPlots.PZ_biomass_stacked, {'Units', 'Position'}, {'inches', [0 0 8 4]})
var = 'PZ_biomass_stacked';

plot_output_timeSeriesPolygon(var, out0, auxVars0, ...
    FixedParams, Forc0, Data0, 'waterMass', waterMass, ...
    'colsP', colsP, 'colsZ', colsZ);


switch savePlots, case true
    % save plots stored in timeSeriesPlots?
    if exist('timeSeriesPlots', 'var')
        fields = fieldnames(timeSeriesPlots);
        for i = 1:length(fields)
            if isvalid(timeSeriesPlots.(fields{i}))
                filename = ['timeSeriesPolygonPlot_' fields{i} '_' waterMass '.png'];
                print(timeSeriesPlots.(fields{i}), fullfile(folder, filename), '-r300', '-dpng');
            end
        end
    end
end

close all
clear timeSeriesPlots


%% Size spectra vs time -- 3D surface plots

waterMass = 'Atlantic';
traj = strcmp(Forc0.waterMass, waterMass);

trophicLevels = {'P','Z'};
nutrientsP = FixedParams.PP_nut;
nutrientsZ = FixedParams.ZP_nut;

Times = 100:FixedParams.nt; % time range to plot
D = FixedParams.nz;

% colourMap = 'plasma';
colourMap = 'jet';

for i = 1:length(trophicLevels)
    trophicLevel = trophicLevels{i};
    N = FixedParams.(['n' trophicLevel 'P_size']);
    Sizes = FixedParams.([trophicLevel 'Pdia']);    
    nutrients = eval(['nutrients' trophicLevel]);
    for j = 1:length(nutrients)
        nutrient = nutrients{j};
        nutrienti = FixedParams.([trophicLevel  'P_' nutrient '_index']);
        y = out0.(trophicLevel);
        y = y(:,:,nutrienti,:,traj);
        y = mean(y, ndims(y), 'omitnan'); % average over trajectories
        y = repmat(reshape(FixedParams.zwidth, [1 D]), [N 1]) .* y; % conc -> quantity
        y = squeeze(sum(y, 2, 'omitnan')); % sum over depths
        y = y(:,Times);

        plotName = [trophicLevel '_' nutrient];
        surfPlots.(plotName) = figure;
        set(surfPlots.(plotName), {'Units', 'Position'}, {'inches', [0 0 7 6]})
        surf(Times, Sizes, y);
        
        set(gca, 'yscale', 'log')
        
        colormap(colourMap);
        shading interp;
        colorbar;
        view(45, 45);
        set(gca, 'XLIM', [min(Times) max(Times)]);
        set(gca, 'YLIM', [min(Sizes) max(Sizes)]);
        set(gcf, 'color', 'white');
        xlabel('year-day')
        ylabel('cell diameter (\mum)')
        switch nutrient
            case 'C', zlab = 'carbon (mmol m^{-2})';
            case 'N', zlab = 'nitrogen (mmol m^{-2})';
            case 'Chl', zlab = 'chlorophyll a (mg m^{-2})';
        end
        zlabel(zlab)
        switch trophicLevel
            case 'P', Title = 'autotrophs';
            case 'Z', Title = 'heterotrophs';
        end
        title(Title)
    end
end



% Surface plots of other quantities integrated or max-value over depth

% Autotroph uptake (production)
RateOrTotal = 'total';

for i = 1:length(nutrientsP)
    nutrient = nutrientsP{i};
    nutrienti = FixedParams.(['PP_' nutrient '_index']);
    
    switch RateOrTotal
        case 'total'
            y = auxVars0.PP_uptake_depthInt;
            y = y(:,:,nutrienti,:,traj);
            y = squeeze(mean(y, ndims(y), 'omitnan')); % average over trajectories
            y = y(:,Times);
            plotName = ['P_uptake_' nutrient];
        case 'rate'
            y = auxVars0.V(1:FixedParams.nPP_size,:,:,:,:);
            y = y(:,:,nutrienti,:,traj);
            y = mean(y, ndims(y), 'omitnan'); % average over trajectories
            y = squeeze(max(y, [], 2, 'omitnan')); % max over depths
            y = y(:,Times);            
            plotName = ['P_uptakeRate_' nutrient];
    end
    
    surfPlots.(plotName) = figure;
    set(surfPlots.(plotName), {'Units', 'Position'}, {'inches', [0 0 7 6]})
    surf(Times, Sizes, y);
    
    set(gca, 'yscale', 'log')
%     set(gca, 'zscale', 'log')
    
    colormap(colourMap);
    shading interp;
    colorbar;
    view(45, 45);
    set(gca, 'XLIM', [min(Times) max(Times)]);
    set(gca, 'YLIM', [min(Sizes) max(Sizes)]);
    set(gcf, 'color', 'white');
    xlabel('year-day')
    ylabel('cell diameter (\mum)')
    switch RateOrTotal
        case 'total'
            switch nutrient
                case 'C'
                    zlab = 'production';
                    Unit = '(mmol C m^{-2} d^{-1})';
                case 'N'
                    zlab = 'uptake';
                    Unit = '(mmol N m^{-2} d^{-1})';
                case 'Chl'
                    zlab = 'production';
                    Unit = '(mg chl a m^{-2} d^{-1})';
            end
        case 'rate'
            Unit = '(d^{-1})';
            switch nutrient
                case {'C', 'Chl'}
                    zlab = 'production rate';
                case 'N'
                    zlab = 'uptake rate';
            end
    end
    zlabel([zlab ' ' Unit])
    title(['autotroph ' zlab])
end


% Predation

% autotroph losses

RateOrTotal = 'total';

for i = 1:length(nutrientsZ)
    nutrient = nutrientsZ{i};
    nutrienti = FixedParams.(['ZP_' nutrient '_index']);
    
    switch RateOrTotal
        case 'total'
            y = auxVars0.predation_losses_depthInt;
            y = y(:,1:FixedParams.nPP_size,:,nutrienti,:,traj);
            y = sum(y); % sum over predators
            y = squeeze(mean(y, ndims(y), 'omitnan')); % average over trajectories
            y = y(:,Times);
            plotName = ['P_grazeLosses_' nutrient];
        case 'rate'
            y = auxVars0.G(:,1:FixedParams.nPP_size,:,:,:);
            y = y(:,:,:,:,traj);
            y = mean(y, ndims(y), 'omitnan'); % average over trajectories
            y = sum(y); % sum over predators
            y = squeeze(max(y, [], 3, 'omitnan')); % max over depths
            y = y(:,Times);            
            plotName = ['P_grazeLossRate_' nutrient];
    end
    
    surfPlots.(plotName) = figure;
    set(surfPlots.(plotName), {'Units', 'Position'}, {'inches', [0 0 7 6]})
    surf(Times, Sizes, y);
    
    set(gca, 'yscale', 'log')
%     set(gca, 'zscale', 'log')
    
    colormap(colourMap);
    shading interp;
    colorbar;
    view(45, 45);
    set(gca, 'XLIM', [min(Times) max(Times)]);
    set(gca, 'YLIM', [min(Sizes) max(Sizes)]);
    set(gcf, 'color', 'white');
    xlabel('year-day')
    ylabel('cell diameter (\mum)')
    switch RateOrTotal
        case 'total'
            zlab = 'grazing losses';
            switch nutrient
                case 'C'
                    Unit = '(mmol C m^{-2} d^{-1})';
                case 'N'
                    Unit = '(mmol N m^{-2} d^{-1})';
            end
        case 'rate'
            Unit = '(d^{-1})';
            zlab = [nutrient ' grazing loss rate'];
    end
    zlabel([zlab ' ' Unit])
    title(['autotroph ' zlab])
end

% heterotroph losses
RateOrTotal = 'total';

for i = 1:length(nutrientsZ)
    nutrient = nutrientsZ{i};
    nutrienti = FixedParams.(['ZP_' nutrient '_index']);
    
    switch RateOrTotal
        case 'total'
            y = auxVars0.predation_losses_depthInt;
            y = y(:,FixedParams.nPP_size+1:end,:,nutrienti,:,traj);
            y = sum(y); % sum over predators
            y = squeeze(mean(y, ndims(y), 'omitnan')); % average over trajectories
            y = y(:,Times);
            plotName = ['Z_grazeLosses_' nutrient];
        case 'rate'
            y = auxVars0.G(:,FixedParams.nPP_size+1:end,:,:,:);
            y = y(:,:,:,:,traj);
            y = mean(y, ndims(y), 'omitnan'); % average over trajectories
            y = sum(y); % sum over predators
            y = squeeze(max(y, [], 3, 'omitnan')); % max over depths
            y = y(:,Times);            
            plotName = ['Z_grazeLossRate_' nutrient];
    end
    
    surfPlots.(plotName) = figure;
    set(surfPlots.(plotName), {'Units', 'Position'}, {'inches', [0 0 7 6]})
    surf(Times, Sizes, y);
    
    set(gca, 'yscale', 'log')
%     set(gca, 'zscale', 'log')
    
    colormap(colourMap);
    shading interp;
    colorbar;
    view(45, 45);
    set(gca, 'XLIM', [min(Times) max(Times)]);
    set(gca, 'YLIM', [min(Sizes) max(Sizes)]);
    set(gcf, 'color', 'white');
    xlabel('year-day')
    ylabel('cell diameter (\mum)')
    switch RateOrTotal
        case 'total'
            zlab = 'grazing losses';
            switch nutrient
                case 'C'
                    Unit = '(mmol C m^{-2} d^{-1})';
                case 'N'
                    Unit = '(mmol N m^{-2} d^{-1})';
            end
        case 'rate'
            Unit = '(d^{-1})';
            zlab = [nutrient ' grazing loss rate'];
    end
    zlabel([zlab ' ' Unit])
    title(['heterotroph ' zlab])
end


% heterotroph gains
RateOrTotal = 'rate';

for i = 1:length(nutrientsZ)
    nutrient = nutrientsZ{i};
    nutrienti = FixedParams.(['ZP_' nutrient '_index']);
    
    y = auxVars0.predation_gains_depthInt;
    y = y(:,:,:,nutrienti,:,traj);
    y = squeeze(mean(y, ndims(y), 'omitnan')); % average over trajectories
    y = y(:,Times);
    plotName = ['Z_grazeGains_' nutrient];
    
    surfPlots.(plotName) = figure;
    set(surfPlots.(plotName), {'Units', 'Position'}, {'inches', [0 0 7 6]})
    surf(Times, Sizes, y);
    
    set(gca, 'yscale', 'log')
%     set(gca, 'zscale', 'log')
    
    colormap(colourMap);
    shading interp;
    colorbar;
    view(45, 45);
    set(gca, 'XLIM', [min(Times) max(Times)]);
    set(gca, 'YLIM', [min(Sizes) max(Sizes)]);
    set(gcf, 'color', 'white');
    xlabel('year-day')
    ylabel('cell diameter (\mum)')
    switch RateOrTotal
        case 'total'
            zlab = 'grazing gains';
            switch nutrient
                case 'C'
                    Unit = '(mmol C m^{-2} d^{-1})';
                case 'N'
                    Unit = '(mmol N m^{-2} d^{-1})';
            end
        case 'rate'
            Unit = '(d^{-1})';
            zlab = [nutrient ' grazing gain rate'];
    end
    zlabel([zlab ' ' Unit])
    title(['heterotroph ' zlab])
end


switch savePlots, case true
    % save plots stored in timeSeriesPlots?
    if exist('surfPlots', 'var')
        fields = fieldnames(surfPlots);
        for i = 1:length(fields)
            if isvalid(surfPlots.(fields{i}))
                filename = ['surfacePlotSizeTime_' fields{i} '_' waterMass '.png'];
                print(surfPlots.(fields{i}), fullfile(folder, filename), '-r300', '-dpng');
            end
        end
    end
end


close all
clear surfPlots


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


switch savePlots, case true
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Older, rougher plotting code below...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Model fit to data

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


switch savePlots, case true
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


switch savePlots, case true
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


switch savePlots, case true
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

