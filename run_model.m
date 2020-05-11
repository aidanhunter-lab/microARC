%% Size-structured 1D monod model with a single predator class

clear
clc
close all

%% Load and filter forcing data
% this is mostly copied from Fabian's code, much of which has been
% commented out

% select your model names
% bioModel  = 'size_based_model'; % What biogeochemistry?
forcModel = 'sinmod'; % What forcing?

% define forcing name (i.e. forcing subfolder)
forcName = 'FramStrait_dep0-100';

% define a run name
runName = 'FramStrait_dep0-100';

% define experiment name
expName = 'AWI-Hausgarten';

% select years
years = 2017:2018;

% set input directories and files
% baseDir = fullfile(filesep, 'media', 'aidan', 'STORAGE', 'Work', 'microARC');
% forcDir = fullfile(baseDir, 'DATA', 'particle_trajectories', forcModel, forcName);
% forcDummy = 'particles_MODEL_DOMAIN_EXPERIMENT-YEAR_YEAR_t*iSub03.mat';

% baseDir = fullfile(filesep, 'media', 'aidan', 'STORAGE', 'Work', 'microARC');
forcDir = fullfile('DATA', 'particle_trajectories', forcModel, forcName);
forcDummy = 'particles_MODEL_DOMAIN_EXPERIMENT-YEAR_YEAR_t*iSub03.mat';


% % set output directory
% outDir = fullfile(baseDir, 'OUTPUT', bioModel, 'FORCSTR', runName, expName);

% select trajectores
%  - selection options
%     1. use all trajectories by setting useTraj = [];
%     2. manually define trajectory index vector, e.g.: iTraj0 = 1:100:7500
%     3. use validation region(s)/transect(s) and time => vector of negative region indices, e.g.: iTraj = -1

% MERGING THIS SCRIPT WITH CODE THAT CHOOSES APPROPRIATE TRAJECTORIES CAN
% BE DONE LATER... FOR NOW JUST SELECT A FEW TRAJECTORIES TO RUN THE MODEL
% OVER

% useTraj = [];
useTraj = 1:100:5000;

% set maximum number of trajectories to be used (set <=0 if no limitation)
nTrajMax = 0;

% list of biological forcing variables (used for initalization)
switch lower(forcModel)
    case 'biomas' % for BIOMAS forcing
        bioForcing = {'NO3', 'Si', 'PS', 'PL', 'ZS', 'ZL', 'ZP'};
    case 'sinmod' % for SINMOD forcing
        bioForcing = {'NO3', 'PS', 'PL'};
end

% outDir = strrep(outDir, 'FORCSTR', sprintf('FRC_%s', forcModel));

% define "validation" region(s) or transect(s); extend list as necessary --
% see Fabian's original code...
% #1: Fram Strait
i = 1;
vali(i).name = 'FramStrait';       % name of region/transect
vali(i).type = 'transect';         % type of validation
vali(i).date = '01-Mar to 30-Sep'; % time when trajectories need to cut the region
vali(i).radius = 25;               % radius around section locations in km

% use this separator string for validation time periods (in the different
% validation regions and/or transects
dateSep = ' to ';

% loop over selected forcing files
for iy = 1:length(years)
    % load forcing
    forcFile = replace(forcDummy, {'MODEL', 'DOMAIN', 'EXPERIMENT', 'YEAR'}, ...
                                  {forcModel, forcName, expName, num2str(years(iy))});
    forcList = dir(fullfile(forcDir, forcFile));
    forcFile = forcList.name;
    try
        tStart = tic;
        load(fullfile(forcDir, forcFile));
        forcing = PDat;
        grid = runDat.grid;
        clear PDat runDat relDat stepsDat
    catch
        error('Forcing file does not exist:\n ==> %s.', fullfile(forcDir, forcFile));
    end
    
    % filter trajectories with geographical jumps
    forcing = filterTrajectoryJumps(forcing);

    % match forcing longitudes with model and mapping coordinate system
    iWest = forcing.x >= 180.0;
    forcing.x(iWest) = forcing.x(iWest) - 360.0;
    clear iWest

    % set model forcing
    nTraj = size(forcing.x,2);
    if  isempty(useTraj)
        % no specific trajectories seleced => use all
        iTraj = 1:nTraj;
    else
        iSel = false(1,nTraj);
        if  sum(useTraj>0)>0
            % used selected trajectory indices
            iSel(useTraj(useTraj>0 & useTraj<=nTraj)) = true;
        end
        iTraj = find(iSel); % convert logical mask into index vector
    end

    % filter forcing data
    y_index = ['y' num2str(years(iy))];
    F.(y_index).iTraj = iTraj;
    F.(y_index).t     = forcing.t;
    F.(y_index).z     = grid.z;
    F.(y_index).x     = forcing.x(:,iTraj);
    F.(y_index).y     = forcing.y(:,iTraj);
    F.(y_index).H     = -forcing.H(:,iTraj);
    F.(y_index).ice   = forcing.ice(:,iTraj);
    F.(y_index).iceh  = forcing.iceh(:,iTraj);
    F.(y_index).swrad = forcing.swrad(:,iTraj);
    F.(y_index).T     = forcing.profiles.temp(:,iTraj,:);
    F.(y_index).kv    = forcing.profiles.Ks(:,iTraj,:) * 0.0001; % cm2/s -> m2/s
    % The 'bioForcing' variables returned from SINMOD can be used to help
    % initialise the model
    for j = 1:length(bioForcing)
        F.(y_index).(sprintf('%sic', bioForcing{j})) = ... 
            forcing.profiles.(bioForcing{j})(:,iTraj,:);
    end
    
    clear forcing
    % Permute trajectory to last dimension
    F.(y_index).T = permute(F.(y_index).T , [1, 3, 2]);
    F.(y_index).kv = permute(F.(y_index).kv , [1, 3, 2]);        
    for j = 1:length(bioForcing)
        F.(y_index).(sprintf('%sic', bioForcing{j})) = ...
            permute(F.(y_index).(sprintf('%sic', bioForcing{j})), [1, 3, 2]);
    end
end



%% Model set-up

% Model input is stored in 3 structs, 'Forc', 'FixedParams' and 'Params',
% that contain the time-varying forcing data, constant parameters, and
% parameters to vary in an optimisation. Any parameter contained in 'Params'
% may be moved into 'FixedParams'.

FixedParams.years = years;
for iy = 1:length(years)
    y_index = ['y' num2str(years(iy))];
    FixedParams.(y_index).nt = length(F.(y_index).t);        % number of timesteps
    FixedParams.(y_index).t = F.(y_index).t;                 % all times
    FixedParams.(y_index).nTraj = length(F.(y_index).iTraj); % number of selected particle trajectories
    FixedParams.(y_index).lat = F.(y_index).y;               % latitude
    FixedParams.(y_index).lon = F.(y_index).x;               % longitude
end


%~~~~~~~~~~~~~~~~~~~~~~
% Modelled depth layers
%~~~~~~~~~~~~~~~~~~~~~~

FixedParams.nz = 15;                                                  % number of modelled depth layers
Htot = 150;                                                           % total modelled depth
dzmin = 5;                                                            % minimum layer width (dzmin <= Htot/(nz-1))
dzmax = 2 * Htot / FixedParams.nz - dzmin;                            % maximum layer width
FixedParams.zw = [0; -cumsum(linspace(dzmin,dzmax,FixedParams.nz))']; % depth of layer edges
FixedParams.zwidth = FixedParams.zw(1:end-1) - FixedParams.zw(2:end); % widths of depth layers
FixedParams.z = 0.5*(FixedParams.zw(1:end-1)+FixedParams.zw(2:end));  % midpoints of depth layers
FixedParams.delz = abs(diff(FixedParams.z));                          % distance between centres of adjacent depth layers


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Prepare forcing data for model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Extracts relavent data and interpolate it over our chosen depth layers
[Forc, FixedParams] = prepareForcing(F,FixedParams);


%~~~~~~~~~~~~~~~~
% State variables
%~~~~~~~~~~~~~~~~

% Inorganic nutrients - only nitrate is modelled
FixedParams.INtype = {'NO3'};
FixedParams.nIN = length(FixedParams.INtype);                               % number of inorganic nutrients

% Plankton
FixedParams.nPP = 6;                                                        % number of phytoplankton size classes
FixedParams.nZP = 1;                                                        % number of zooplankton classes

% Phytoplanton sizes - smallest diameter is 0.5 mu m, volumes of successive size 
% classes increase by factors of 32 (equally spaced on log-scale).
PPdia = 0.5;
FixedParams.PPsize = 4/3*pi*(PPdia/2)^3;
FixedParams.PPsize(2:FixedParams.nPP) = 32 .^ (1:FixedParams.nPP-1) .* ...
    FixedParams.PPsize(1);                                                  % cell volumes [mu m^3]
PPdia = 2 .* (3 .* FixedParams.PPsize ./ (4*pi)) .^ (1/3);
FixedParams.diatoms = FixedParams.PPsize >= 100;                            % assume large phytoplankton are diatoms - only needed to split SINMOD output over classes during state variable initialisation

FixedParams.phytoplankton = [true(1,FixedParams.nPP) ... 
    false(1,FixedParams.nZP)];                                              % index phytoplankton
FixedParams.zooplankton = [false(1,FixedParams.nPP) ... 
    true(1,FixedParams.nZP)];                                               % index zooplankton

% Organic matter - only DOM is explicitly modelled
FixedParams.OMtype = {'DOM'};
FixedParams.nOM = length(FixedParams.OMtype);                               % number of organic matter variables

FixedParams.nVar = FixedParams.nIN + FixedParams.nPP + FixedParams.nZP + ...
    FixedParams.nOM;                                                        % number of state variables per depth and trajectory


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% State variable initial condition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Store state variables in array: 1st dimension = variable
%                                 2nd dimension = location (trajectory)
% 
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth]
%                      zooplankton         [depth]
%                      organic matter      [depth]

[v0, initialState] = initialiseVariables(FixedParams,Forc); % v0 stores initial input vectors; initialState contains the same info in a more readable struct


%~~~~~~~~~~~~~~~~~
% Variable indexes
%~~~~~~~~~~~~~~~~~

% Create indexes to extract inorganic nutrients, plankton and organic
% matter from the array containing all variables, and indexes for specific
% nutrient types within each of IN, PP, ZP and OM.

FixedParams = createIndexes(FixedParams);



%~~~~~~~~~~~~~~~~~
% Model parameters
%~~~~~~~~~~~~~~~~~

FixedParams.POM_is_lost = true; % is POM lost from the system by sinking below bottom modelled depth layer

% Initialise model parameters
[Params, FixedParams] = initialiseParameters(FixedParams);
Params0 = Params;



%% Integration time interval
dt=1;             % output timestep (days)
t0=1;             % tstep 1
tspan = nan(length(FixedParams.years),2);
tmax = nan(1,length(FixedParams.years));
for iy = 1:length(FixedParams.years)
    y_index = ['y' num2str(FixedParams.years(iy))];
    tmax(iy) = FixedParams.(y_index).nt;  % final timestep
    tspan(iy,:) = [t0 tmax(iy)];          % modelled time interval
end


%% Integrate

nEquations = FixedParams.nVar * FixedParams.nz; % number of ODEs

% solve with ode45 solver
odeDTmin = 1/24/60; % min timestep = 1 minute
odeDTmax = 1;      % max timestep = 1 day

% set solver options: positive definite, tolerances, min/max tsteps
ode45options=odeset('NonNegative',[1 ones(1, nEquations)],...
  'AbsTol',1e-6,'RelTol',1e-4,...
  'InitialStep',odeDTmin,'MaxStep',odeDTmax);


% Storage for output
for iy = 1:length(FixedParams.years)
    y_index = ['y' num2str(FixedParams.years(iy))];
    out.(y_index) = nan(nEquations, tmax(iy), FixedParams.(y_index).nTraj);
end

% Optional extra outputs to extract - 4 cases allowed so far. When running
% a parameter optimisation algorithm set returnExtras='none' for speed.
returnExtras = 'none';
% returnExtras = 'auxiliaryAndRates';
%returnExtras = 'auxiliary';             % extra, non-state variables, are stored
% returnExtras = 'rates';               % state variable rates of change are stored

FixedParams.returnExtras = strcmp(returnExtras, 'auxiliary') || ...
    strcmp(returnExtras, 'auxiliaryAndRates');

% create storage for extra outputs if required
switch returnExtras
    case 'auxiliaryAndRates'
        for iy = 1:length(FixedParams.years)
            y_index = ['y' num2str(FixedParams.years(iy))];
            auxVars.(y_index) = [];
            rates.(y_index) = nan(nEquations, tmax(iy), FixedParams.(y_index).nTraj);
        end
    case 'auxiliary'
        for iy = 1:length(FixedParams.years)
            y_index = ['y' num2str(FixedParams.years(iy))];
            auxVars.(y_index) = [];
        end
    case 'rates'
        for iy = 1:length(FixedParams.years)
            y_index = ['y' num2str(FixedParams.years(iy))];
            rates.(y_index) = nan(nEquations, tmax(iy), FixedParams.(y_index).nTraj);
        end
end


% Combine parameters into single input structure
parameterList.Params = Params;
parameterList.FixedParams = FixedParams;


% Loop through years
for iy = 1:length(FixedParams.years)    
    y_index = ['y' num2str(FixedParams.years(iy))];    
    % Loop through trajectories
    for kk = 1:FixedParams.(y_index).nTraj
        % Extract forcing data for year iy, trajectory kk; and permute time to
        % last dimension
        forcing.T = permute(Forc.(y_index).T(:,:,kk), [2 1]);
        forcing.PAR = permute(Forc.(y_index).PAR(:,:,kk), [2 1]);
        forcing.K = permute(Forc.(y_index).kv(:,:,kk), [2 1]);
        
        % Initial condition for year iy, trajectory kk
        v_in = v0.(y_index)(:,kk);        
        out.(y_index)(:,1,kk) = v_in;
        
        % Initialise extra outputs
        if ~strcmp(returnExtras, 'none')
            ts = tspan(iy,1):tspan(iy,1)+1;
            parameterList.Forc.T = forcing.T(:,ts);
            parameterList.Forc.PAR = forcing.PAR(:,ts);
            parameterList.Forc.K = forcing.K(:,ts);
            switch returnExtras
                case 'auxiliaryAndRates'
                    [dvdt,extraOutput] = ODEs(0, v_in, parameterList);
                    fields = fieldnames(extraOutput);
                    for jj = 1:length(fields)
                        auxVars.(y_index).(fields{jj})(:,1,kk) = extraOutput.(fields{jj});
                    end
                    rates.(y_index)(:,1,kk) = dvdt;
                case 'auxiliary'
                    [~,extraOutput] = ODEs(0, v_in, parameterList);
                    fields = fieldnames(extraOutput);
                    for jj = 1:length(fields)
                        auxVars.(y_index).(fields{jj})(:,1,kk) = extraOutput.(fields{jj});
                    end
                case 'rates'
                    dvdt = ODEs(0, v_in, parameterList);
                    rates.(y_index)(:,1,kk) = dvdt;
            end
            clear extraOutput dvdt
        end
        
        % Integrate step-wise between successive data points -- perhaps the
        % only way to guarentee smoothly continuous ODE functions
        tic        
%         profile on
        for tt = 1:tmax(iy)-1            
            % Set forcing data
            parameterList.Forc.T = forcing.T(:,tt:tt+1);
            parameterList.Forc.PAR = forcing.PAR(:,tt:tt+1);
            parameterList.Forc.K = forcing.K(:,tt:tt+1);
            
            % Integrate -- returning full output structure can be useful for debugging
            sol=ode45(@(t, v_in) ODEs(t, v_in, parameterList), [0 1], v_in, ode45options);
            
            % Store solutions each day (each forcing data time-step)
            out.(y_index)(:,tt+1,kk) = deval(sol, 1);            
            clear sol
            
            % Update initials for next time step
            v_in = out.(y_index)(:,tt+1,kk);
            
            % Extract extra outputs            
            if ~strcmp(returnExtras, 'none')
                switch returnExtras
                    case 'auxiliaryAndRates'
                        [dvdt,extraOutput] = ODEs(0, v_in, parameterList);
                        fields = fieldnames(extraOutput);
                        for jj = 1:length(fields)
                            auxVars.(y_index).(fields{jj})(:,tt+1,kk) = extraOutput.(fields{jj});
                        end
                        rates.(y_index)(:,1,kk) = dvdt;
                    case 'auxiliary'
                        [~,extraOutput] = ODEs(0, v_in, parameterList);
                        fields = fieldnames(extraOutput);
                        for jj = 1:length(fields)
                            auxVars.(y_index).(fields{jj})(:,tt+1,kk) = extraOutput.(fields{jj});
                        end
                    case 'rates'
                        dvdt = ODEs(0, v_in, parameterList);
                        rates.(y_index)(:,tt+1,kk) = dvdt;
                end
                clear extraOutput dvdt
            end
        end
%         profile viewer
        toc
        
    end    
end



%% Extract outputs

% State variables
for iy = 1:length(FixedParams.years)
    y_index = ['y' num2str(FixedParams.years(iy))];
    output.(y_index).N = out.(y_index)(FixedParams.IN_index,:,:);    
    output.(y_index).P = reshape(out.(y_index)(FixedParams.PP_index,:,:), ...
        [FixedParams.nPP FixedParams.nz FixedParams.(y_index).nt ...
        FixedParams.(y_index).nTraj]);
    output.(y_index).Z = out.(y_index)(FixedParams.ZP_index,:,:);
    output.(y_index).OM = out.(y_index)(FixedParams.OM_index,:,:);
end


%% Check mass is conserved

% This will need adjusted for case POM_is_lost==true
for iy = 1:length(FixedParams.years)
    if iy == 1, figure, end
    y_index = ['y' num2str(FixedParams.years(iy))];    
    totalMass.(y_index) = squeeze(sum(FixedParams.zwidth .* output.(y_index).N)) + ...
        squeeze(sum(FixedParams.zwidth .* output.(y_index).Z)) + ...
        squeeze(sum(FixedParams.zwidth .* output.(y_index).OM)) + ...    
        squeeze(sum(sum(FixedParams.zwidth .* permute(output.(y_index).P, [2 1 3 4])))); % total quantity at each time step - assuming 1m^2 cross-sectional area of water column
    subplot(length(FixedParams.years), 1, iy)
    plot(diff(totalMass.(y_index)))
    if iy == length(FixedParams.years)
        xlabel('time step')
        suptitle('change in total mass each time step')
    end
    ylabel('mass difference')
    title(num2str(FixedParams.years(iy)))
end




%% Plot some output...

close all

iy = 1;    % Choose year
y_index = ['y' num2str(FixedParams.years(iy))];
nsimTraj = size(out.(y_index), 3);
nt = size(out.(y_index), 2);

kk = 1;  % Choose output from a single trajectory

N = output.(y_index).N(:,:,kk);
P = output.(y_index).P(:,:,:,kk);
Z = output.(y_index).Z(:,:,kk);
OM = output.(y_index).OM(:,:,kk);

% Save plots?
save = false;
folder = 'OUTPUT/plots/';

% Create time-depth grid for interpolation
[depth, time] = ndgrid(abs(FixedParams.z), 1:FixedParams.(y_index).nt);
[depthGrid, timeGrid] = ndgrid(1:1:abs(FixedParams.zw(end)), ...
    1:FixedParams.(y_index).nt);

% Should plots be smoothed by linear interpolation?
smooth = 'nearest'; % interpolation type
% smooth = 'linear';

% Plot forcing data
figure(1)

subplot(3,1,1)
x = Forc.(y_index).T(:,:,kk)';
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
contourf(Fsmooth)
colorbar
title('Temperature (\circC)')
ylabel('depth (m)')
xticks(100:100:FixedParams.(y_index).nt)
xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
yticks(linspace(0,abs(FixedParams.zw(end)),7))
yticklabels(linspace(FixedParams.zw(end),0,7))

subplot(3,1,2)
x = Forc.(y_index).kv_center(:,:,kk)';
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
contourf(Fsmooth)
colorbar
title('Diffusivity (m^2 day^{-1})')
ylabel('depth (m)')
xticks(100:100:FixedParams.(y_index).nt)
xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
yticks(linspace(0,abs(FixedParams.zw(end)),7))
yticklabels(linspace(FixedParams.zw(end),0,7))

subplot(3,1,3)
x = Forc.(y_index).PAR(:,:,kk)';
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
contourf(Fsmooth)
colorbar
title('PAR (\muEin day^{-1} m^{-2})')
xlabel('year-day')
ylabel('depth (m)')
xticks(100:100:FixedParams.(y_index).nt)
xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
yticks(linspace(0,abs(FixedParams.zw(end)),7))
yticklabels(linspace(FixedParams.zw(end),0,7))

colormap plasma

if save
    filename = 'forcing_data.png';
    figFile = fullfile(folder, filename);    
    fig = gcf;
    % Adjust figure window dimensions
    fig.Units = 'inches';
    fig.Position = [0 0 6 7];
    fig.PaperPositionMode = 'auto'; % save figure window as-is
    print(fig, figFile, '-r300', '-dpng');
end



% Inorganic nutrient
figure(2)
x = N;
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
contourf(Fsmooth)
colorbar
title('DIN (mmol N / m^3)')
xlabel('year-day')
ylabel('depth (m)')
xticks(100:100:FixedParams.(y_index).nt)
xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
yticks(linspace(0,abs(FixedParams.zw(end)),7))
yticklabels(linspace(FixedParams.zw(end),0,7))
colormap plasma

if save
    filename = 'dissolved_inorganic_nutrient.png';
    figFile = fullfile(folder, filename);    
    fig = gcf;
    % Adjust figure window dimensions
    fig.Units = 'inches';
    fig.Position = [0 0 8 4];
    fig.PaperPositionMode = 'auto'; % save figure window as-is
    print(fig, figFile, '-r300', '-dpng');
end



figure(3)
multiPanelPlot = strcmp(returnExtras, 'auxiliary') || ...
    strcmp(returnExtras, 'auxiliaryAndRates');
if multiPanelPlot, subplot(3,1,1), end
x = OM;
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
contourf(Fsmooth)
colorbar
title('DON (mmol N / m^3)')
if ~multiPanelPlot, xlabel('year-day'), end
ylabel('depth (m)')
xticks(100:100:FixedParams.(y_index).nt)
xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
yticks(linspace(0,abs(FixedParams.zw(end)),7))
yticklabels(linspace(FixedParams.zw(end),0,7))

if multiPanelPlot
    subplot(3,1,2)
    
    x = auxVars.(y_index).POM(:,:,kk);
    F = griddedInterpolant(depth, time, x, smooth);
    Fsmooth = flip(F(depthGrid, timeGrid));
    contourf(Fsmooth)
    colorbar
    title('PON before sinking (mmol N / m^3)')
    ylabel('depth (m)')
    xticks(100:100:FixedParams.(y_index).nt)
    xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
    yticks(linspace(0,abs(FixedParams.zw(end)),7))
    yticklabels(linspace(FixedParams.zw(end),0,7))
    
    subplot(3,1,3)
    x = auxVars.(y_index).remin_POM(:,:,kk);
    F = griddedInterpolant(depth, time, x, smooth);
    Fsmooth = flip(F(depthGrid, timeGrid));
    contourf(Fsmooth)
    colorbar
    title('PON after sinking (mmol N / m^3)')
    xlabel('year-day')
    ylabel('depth (m)')
    xticks(100:100:FixedParams.(y_index).nt)
    xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
    yticks(linspace(0,abs(FixedParams.zw(end)),7))
    yticklabels(linspace(FixedParams.zw(end),0,7))
    
end

colormap plasma

if save
    filename = 'organic_nutrient.png';
    figFile = fullfile(folder, filename);    
    fig = gcf;
    % Adjust figure window dimensions
    fig.Units = 'inches';
    if ~multiPanelPlot
        fig.Position = [0 0 5 4];
    else
        fig.Position = [0 0 6 7];        
    end
    fig.PaperPositionMode = 'auto'; % save figure window as-is
    print(fig, figFile, '-r300', '-dpng');
end






figure(4) % phytoplankton - different panel for each size

nr = floor(sqrt(FixedParams.nPP));
nc = ceil(FixedParams.nPP / nr);

for ii = 1:FixedParams.nPP
    subplot(nr,nc,ii)
    x = squeeze(P(ii,:,:));
    F = griddedInterpolant(depth, time, x, smooth);
    Fsmooth = flip(F(depthGrid, timeGrid));
    contourf(Fsmooth)
    colorbar
    title([num2str(round(FixedParams.PPsize(ii),2,'significant')) ' \mum^3'])
    xlabel('year-day')
    ylabel('depth (m)')
    xticks(100:100:FixedParams.(y_index).nt)
    xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
    yticks(linspace(0,abs(FixedParams.zw(end)),7))
    yticklabels(linspace(FixedParams.zw(end),0,7))
end
suptitle('phytoplankton abundance (mmol N / m^3)')
colormap plasma

if save
    filename = 'phytoplankton.png';
    figFile = fullfile(folder, filename);    
    fig = gcf;
    % Adjust figure window dimensions
    fig.Units = 'inches';
    fig.Position = [0 0 12 6];
    fig.PaperPositionMode = 'auto'; % save figure window as-is
    print(fig, figFile, '-r300', '-dpng');
end



figure(5) % zooplankton

x = Z;
F = griddedInterpolant(depth, time, x, smooth);
Fsmooth = flip(F(depthGrid, timeGrid));
contourf(Fsmooth)
colorbar
title('zooplankton abundance (mmol N / m^3)')
xlabel('year-day')
ylabel('depth (m)')
xticks(100:100:FixedParams.(y_index).nt)
xticklabels(yearday(FixedParams.(y_index).t(100:100:FixedParams.(y_index).nt)))
yticks(linspace(0,abs(FixedParams.zw(end)),7))
yticklabels(linspace(FixedParams.zw(end),0,7))
colormap plasma

if save
    filename = 'zooplankton.png';
    figFile = fullfile(folder, filename);    
    fig = gcf;
    % Adjust figure window dimensions
    fig.Units = 'inches';
    fig.Position = [0 0 5 4];
    fig.PaperPositionMode = 'auto'; % save figure window as-is
    print(fig, figFile, '-r300', '-dpng');
end





