%% Size-structured 1D NPZD model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% run model on a parameter ensemble
% in both, summer and autumn configuration

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Refresh the workspace
clear; clc; close all; delete(gcp('nocreate'))

% add path of LHS sets
addpath(genpath('../Documents/microARC model/'))
% Include microARC path
addpath(genpath(fileparts(which('run_model'))))


% Store folders/filenames of data and saved parameters
mkdir('results/VanessasEnsemble/seasonalConfig/')
mkdir('results/VanessasEnsemble/plots/')

Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', []);
display(Directories)
Directories.resultsDir  = 'results/VanessasEnsemble/seasonalConfig/'
Directories.plotDir  = 'results/VanessasEnsemble/plots/'

% set up model
% Use default set-up to load and organise data and initialise model parameters.
% Set-up may be modified here by passing some extra arguments as name-value
% pairs (preferable), or directly modified within modelSetUp.m
% summer
[Forc_s, FixedParams_s, Params_s, Data_s] = modelSetUp2(Directories, ...
    'displayAllOutputs', true, 'seasonConfig', 'summer', 'years', 2018) 
% autumn
[Forc_a, FixedParams_a, Params_a, Data_a] = modelSetUp2(Directories, ...
    'displayAllOutputs', true, 'seasonConfig', 'autumn', 'years', 2018, ...
    'maxDist',100) % maxDist is 100 here, because when using default (25) all trajs would be of atlantic origin

% Params are equal for summer and autumnConfig. 
Params = Params_s;
clear Params_s Params_a 

% Useful name-value pairs for modelSetUp.m include: useTraj, numTraj, ESDmin, ESDmax, nsizes

% modelSetup.m may also produce plots -- set to 'true' any name-value pair as shown below

% [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
%     'plotCellConcSpectra', true, 'plotBioVolSpectra', true, ...
%     'plotSizeClassIntervals', true, ...
%     'trajectoryPlot', true, 'dendrogramPlot', true, ...
%     'plotScalarData', true, 'plotSizeData', true, 'plotAllData', true, ...
%     'displayForc', true, 'displayData', true);


% map the tjaectories and choose one arctic and one atlantic
map_trajectories(Forc_s, FixedParams_s);
% id 2935 3126 for arc, 621 883 atl
map_trajectories(Forc_a, FixedParams_a);
% id 291 298 for arc, 220 246 atl


% OPTIONAL: filter out single trajectories for shorter model run time
singleTrajSetUp = true;

if singleTrajSetUp
    
    Forc_s_ = Forc_s;
    fieldnames = fields(Forc_s);

    % summerTrajs = [621 883 2935 3126];
    summerTrajs = [343 3126];
    
    % summerTrajsIndex = [find(Forc_s.iTraj == summerTrajs(1)) find(Forc_s.iTraj == summerTrajs(2))];
    summerTrajsIndex = find(ismember(Forc_s.iTraj, summerTrajs));  % also need later for obersvation-modOutput matching
    nTrajOrig = Forc_s.nTraj;

    for i = 1:length(fieldnames)

        if size(Forc_s.(fieldnames{i}),2) == nTrajOrig
            Forc_s_.(fieldnames{i}) = Forc_s.(fieldnames{i})(:,summerTrajsIndex);
        elseif size(Forc_s.(fieldnames{i}),3) == nTrajOrig
            Forc_s_.(fieldnames{i}) = Forc_s.(fieldnames{i})(:,:,summerTrajsIndex);
        end
    end
    Forc_s_.nTraj = length(summerTrajs);
    Forc_s = Forc_s_;
    
    %autumn
    Forc_a_ = Forc_a;
    fieldnames = fields(Forc_a);

    % autumnTrajs = [220 246 291 298];
    autumnTrajs = [220 294]; 
    % autumnTrajsIndex = [find(Forc_a.iTraj == autumnTrajs(1)) find(Forc_a.iTraj == autumnTrajs(2))];
    autumnTrajsIndex = find(ismember(Forc_a.iTraj, autumnTrajs));
    nTrajOrig = Forc_a.nTraj;

    for i = 1:length(fieldnames)

        if size(Forc_a.(fieldnames{i}),2) == nTrajOrig
            Forc_a_.(fieldnames{i}) = Forc_a.(fieldnames{i})(:,autumnTrajsIndex);
        elseif size(Forc_a.(fieldnames{i}),3) == nTrajOrig
            Forc_a_.(fieldnames{i}) = Forc_a.(fieldnames{i})(:,:,autumnTrajsIndex);
        end
    end
    Forc_a_.nTraj = length(autumnTrajs);
    Forc_a = Forc_a_;

    clear Forc_s_ Forc_a_ fieldnames nTrajOrig
    
    % map again to check if selection was successful
    map_trajectories(Forc_s, FixedParams_s);
    % id 2935 for arc, 621 atl
    map_trajectories(Forc_a, FixedParams_a);
    % id 291 for arc, 220 atl
end


% params to be tuned:
% A, h, m2, aP, theta, xi, aG, sigG, Lambda, lambda_max, rDON (=rDOC), rPOM
% (= POC), 
% Qmin_QC_a, Qmin_QC_b, Qmax_delQ_a, Qmax_delQ_b, Vmax_QC_a, Vmax_QC_b
% aN_QC_a, aN_QC_b, pmax_a, pmax_b, Gmax_a, Gmax_b, m_a, m_b, wp_a, wp_b,
% beta1, beta2, beta3, wPOM1
% reduce this number by excluding parameters that show to be insensitive
% after first analysis




% get names and bounds of all non-fixed parameters
Pnames = { 'A', 'h', 'm2', 'aP', 'theta', 'xi', 'aG', 'sigG', 'Lambda', ...
    'lambda_max','Qmin_QC_a', 'Qmin_QC_b', 'Qmax_delQ_a',...
    'Qmax_delQ_b', 'Vmax_QC_a', 'Vmax_QC_b', 'aN_QC_a', 'aN_QC_b',...
    'pmax_a', 'pmax_b', 'Gmax_a', 'Gmax_b', 'm_a', 'm_b', 'wp_a', 'wp_b',...
    'beta1', 'beta2', 'beta3', 'wPOM1', 'rDON' , 'rPON', 'rDOC', 'rPOC'};
% get bounds of all tunable Params
for i = 1:length(Pnames)
    
    Pmin(i) = Params.bounds.(Pnames{i})(1);
    Pmax(i) = Params.bounds.(Pnames{i})(2); 
end

% rDOC and rPOC equal rDON and rPON, so add them at the end. Don't forget
% the respective columns in the LHC matrix! 

% read LHC matrix
LHS = load([Directories.resultsDir 'LHS_p32_2000.txt']) ;

ncol = size(LHS,2);
% clone last two colomns for rDOC and rPOC
LHS(:, [ncol+1 ncol+2]) = LHS(:, [ncol-1 ncol]); 

deltaP = Pmax - Pmin;

Pensemble = Pmin + LHS .* deltaP;  % each row represents the parameter set for one ensemble member (2000 sets of 32+2 params)

% save Pensemble:
save [Directories.resultsDir 'parameterEnsemble2000.txt'] Pensemble -ascii


% preparation of ensemble parameter sets is now complete. 

%% Run ensemble model
% open a loop, repeat for each ensemble member (nrow(Pensemble))
% ... set up model, directories, all that

% second: updateParameters, run the model, and
% save the results as modData to save memory. 
% each iteration (member) is using one parameter set 
% from the parameter ensemble


errorMembers_s = [];
errorMembers_a = [];

tStart = tic; 
for i = 1:size(Pensemble,1)

    tic; disp(['.. started ensemble member ' num2str(i)  ' at']); disp(datetime('now'))
    
    
    % important step: update parameters 
    % Parameters may be modified using name-value pairs in updateParameters.m
    % Params = updateParameters(Params, FixedParams, 'pmax_a', 20, 'aP', 0.01);
    
    % create namelist of parameter name and value
    ensArgIn = {};
    for j = 1:length(Pnames)
       ensArgIn{1,j} = Pnames{j};
       ensArgIn{2,j} = Pensemble(i,j);  % i with index of this ensemble member, j of parameter
    end
    
    % pass namelist to updateParameters function
    ensParams = updateParameters(Params, FixedParams_s, ensArgIn{:});
    

    
    % run model
    % Input data (forcing trajectories and fitting data) may be filtered by the
    % origin of particle trajectories (Atlantic or Arctic)
    % [Forc, Data] = filterInputByOrigin(Forc, Data, 'fitTrajectories', 'Atlantic');

    % Integrate
    % Initialise state variables.
    % Store state variables in array v0: 1st dimension = variable
    %                                    2nd dimension = location (trajectory)
    v0_s = initialiseVariables(FixedParams_s, ensParams, Forc_s);
    v0_a = initialiseVariables(FixedParams_a, ensParams, Forc_a);
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
    % tic; disp('.. started integration at'); disp(datetime('now'))
    
    try
        [out_s, auxVars_s] = integrateTrajectories(FixedParams_s, ensParams, Forc_s, v0_s, ... 
         FixedParams_s.odeIntegrator, FixedParams_s.odeOptions);
    catch
        warning(['Error in integration of ensemble member ' i, ... 
            'skipping to next iteration'])
        errorMembers_s = [errorMembers_s, i] 
        continue
    end
    
    try 
    [out_a, auxVars_a] = integrateTrajectories(FixedParams_a, ensParams, Forc_a, v0_a, ... 
        FixedParams_a.odeIntegrator, FixedParams_a.odeOptions);
    catch
        warning(['Error in integration of ensemble member ' i, ... 
            'skipping to next iteration'])
        errorMembers_a = [errorMembers_a, i] 
        continue
    end
 
    % display(out)
    % display(auxVars)
    
    % save output (name should corresponf to ensemble ID)
    % saving out creates a very large output
    % save out.m out % produces 255 mb; with 2000 runs 500 GB of storage would be needed
    % we need to condense the output
    % 1. using the modData format: only saves model output where in-situ
    % data are available; for cost calculation

    modData_s = matchModOutput2Data2(out_s, auxVars_s, Data_s, FixedParams_s, ...
        'fitToFullSizeSpectra', true, 'Forc', Forc_s, 'evTrajIX', summerTrajsIndex);
    modData_a = matchModOutput2Data2(out_a, auxVars_a, Data_a, FixedParams_a, ...
        'fitToFullSizeSpectra', true, 'Forc', Forc_a, 'evTrajIX', autumnTrajsIndex);
    % cellDensity and BiovolumeDensity are now extracted for full size 
    % spectra of 10 trajectories per sampling event. 

    
    
    % 2. average over (all?) trajectories...
    % this is for he validation against the WOA and remote sensing data
    % extract4ed aling the model trajectories
    
    % create meanOUT (over traj dimension)
    % summer 
    meanOut_s.N.atlantic = mean(out_s.N(:,:,:,strcmp(Forc_s.waterMass, 'Atlantic')),4);
    meanOut_s.N.arctic = mean(out_s.N(:,:,:,strcmp(Forc_s.waterMass, 'Arctic')),4);
    
    meanOut_s.P.atlantic = mean(out_s.P(:,:,:,:,strcmp(Forc_s.waterMass, 'Atlantic')),5);
    meanOut_s.P.arctic = mean(out_s.P(:,:,:,:,strcmp(Forc_s.waterMass, 'Arctic')),5);
    
    meanOut_s.Z.atlantic = mean(out_s.Z(:,:,:,:,strcmp(Forc_s.waterMass, 'Atlantic')),5);
    meanOut_s.Z.arctic = mean(out_s.Z(:,:,:,:,strcmp(Forc_s.waterMass, 'Arctic')),5);
    
    meanOut_s.OM.atlantic = mean(out_s.OM(:,:,:,:,strcmp(Forc_s.waterMass, 'Atlantic')),5);
    meanOut_s.OM.arctic = mean(out_s.OM(:,:,:,:,strcmp(Forc_s.waterMass, 'Arctic')),5);
    % autumn 
    meanOut_a.N.atlantic = mean(out_a.N(:,:,:,strcmp(Forc_a.waterMass, 'Atlantic')),4);
    meanOut_a.N.arctic = mean(out_a.N(:,:,:,strcmp(Forc_a.waterMass, 'Arctic')),4);
    
    meanOut_a.P.atlantic = mean(out_a.P(:,:,:,:,strcmp(Forc_a.waterMass, 'Atlantic')),5);
    meanOut_a.P.arctic = mean(out_a.P(:,:,:,:,strcmp(Forc_a.waterMass, 'Arctic')),5);
    
    meanOut_a.Z.atlantic = mean(out_a.Z(:,:,:,:,strcmp(Forc_a.waterMass, 'Atlantic')),5);
    meanOut_a.Z.arctic = mean(out_a.Z(:,:,:,:,strcmp(Forc_a.waterMass, 'Arctic')),5);
    
    meanOut_a.OM.atlantic = mean(out_a.OM(:,:,:,:,strcmp(Forc_a.waterMass, 'Atlantic')),5);
    meanOut_a.OM.arctic = mean(out_a.OM(:,:,:,:,strcmp(Forc_a.waterMass, 'Arctic')),5);

    
    % save output 
    ensMemberName = sprintf( '%04d', i );
    
    % summer
    filename = [Directories.resultsDir ensMemberName '_EnsembleMemberOutput_summer.m'];
    m = matfile(filename, 'Writable', true);
    m.modData = modData_s;
    m.meanOut = meanOut_s; 
    m.fixedParams = FixedParams_s;
    m.ensParams = ensParams; 
    % clear output
     clear out_s, auxVars_s, modData_s, meanOut_s
    % autumn
    filename = [Directories.resultsDir ensMemberName '_EnsembleMemberOutput_autumn.m'];
    m = matfile(filename, 'Writable', true);
    m.modData = modData_a;
    m.meanOut = meanOut_a; 
    m.fixedParams = FixedParams_a;
    m.ensParams = ensParams; 
    % clear output
     clear out_a, auxVars_a, modData_a, meanOut_a
     

    clear ensParams
    
    disp(['.. ended ensemble member ' num2str(i)  ' at']); disp(datetime('now')); toc
    
end
disp(['.. ended ensemble run ' num2str(i)  ' at']); disp(datetime('now')); toc(tStart)


% then continue with model cost calculation in
% getEnsembleCost_summerAutumn.m


