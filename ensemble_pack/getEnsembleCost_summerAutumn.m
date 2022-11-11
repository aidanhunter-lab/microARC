%% Calculate fit of model ensemble to observational data


% Refresh the workspace
clear; clc; close all; delete(gcp('nocreate'))

% add path of LHS sets
addpath(genpath('../Documents/microARC model/'))
% Include microARC path
addpath(genpath(fileparts(which('run_model'))))


Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', []);
Directories.resultsDir  = 'results/VanessasEnsemble/seasonalConfig/';
Directories.plotDir  = 'results/VanessasEnsemble/plots/';
display(Directories)

% load (default) model set up, including Data
[Forc_s, FixedParams_s, Params_s, Data_s] = modelSetUp2(Directories, ...
    'displayAllOutputs', true, 'seasonConfig', 'summer', 'years', 2018) 
% autumn
[Forc_a, FixedParams_a, Params_a, Data_a] = modelSetUp2(Directories, ...
    'displayAllOutputs', true, 'seasonConfig', 'autumn', 'years', 2018, ...
    'maxDist',100) % maxDist is 100 here, because when using default (25) all trajs would be of atlantic origin




% Params are equal for summer and autumnConfig. 
Params = Params_s;
clear Params_s Params_a 

% List available cost function choices
[~,~,~,cfc] = costFunction2();
% Select cost function from cfc
disp(cfc)



%% begin cost calculation 
% skip and load results if no changes were done

ensembleCost = [];  % make it a struct instead of a cell!

seasonConfigs = {'summer', 'autumn'}

for i= 1:2000 % for each enselmble member
   
    for is = 1:length(seasonConfigs)
        sc = seasonConfigs{is}; 
      
        % load ensemble output (modData)
        ensOut = matfile([Directories.resultsDir, sprintf( '%04d', i ) ,'_EnsembleMemberOutput_',sc , '.m']);


        % select correct Data for model-data-comparision
        switch sc
            case 'summer' 
                Data = Data_s;
            case 'autumn' 
                Data = Data_a; 
        end

%         % calculate cost of each ensemble member

        % IQD
        % based on station spectra
        [costIQD, costComponentsIQD] = costFunction2('label', 'IQD_vectorFullSpectrum', ...         
            'Data', Data, 'modData', ensOut.modData); 

        % based on scenario spetra
        % was able to convince Markus to abandon this approach :)
%         [costIQD_S, costComponentsIQD_S] = costFunction2('label', 'IQD_S_vectorFullSpectrum', ...         
%             'Data', Data, 'modData', ensOut.modData);                                                       
%         
        
        % Synthetic Likelihood
        [costL, costComponentsL] = costFunction2('label', 'syntheticLikelihood_scalar', ...
            'Data', Data, 'modData', ensOut.modData); 
        
%         % Hellinger for FullSpectrum, no grouping (for comparision
%         against IQD)
%         [costHD, costComponentsHD] = costFunction2('label', 'Hellinger_vectorFullSpectrum', ...
%             'Data', Data, 'modData', ensOut.modData);
  


        ensembleCost.(sc).IQD(i) = costIQD;  % restructure ensembleCost as struct instead of cell array
%        ensembleCost.(sc).IQD_S(i) = costIQD_S;
        ensembleCost.(sc).L(i) = costL;
%        ensembleCost.(sc).HD(i) = costHD;
     
        % display(i)
    end
end




% save ensembleCost
%filename = [Directories.resultsDir 'ensembleCost_binned.m'];
filename = [Directories.resultsDir 'ensembleCost_summerAutumn_lightlim_diffu.m'];
%     m = matfile(filename, 'Writable', true);
%     m.cost = ensembleCost;


% load ensembleCost
costfile = matfile(filename)
ensembleCost = costfile.cost; 



%% some plotting
% not beautiful at all, just for exploring the output
% some code is very old and likeli not functional anymore

% preview plot

% summer
plt = figure
s = scatter(ensembleCost.summer.IQD, ensembleCost.summer.L)
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (IQD\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for all ensemble members, summer config')


% autumn
plt = figure
s = scatter(ensembleCost.autumn.IQD, ensembleCost.autumn.L)
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (IQD\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for all ensemble members, autumn config')


% 
% % erste Scatterplots:
% plt = figure
% s = scatter(ensembleCost(1,:), ensembleCost(2,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:1475);
% 
% plt = figure
% s = scatter(ensembleCost(3,:), ensembleCost(4,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:1475);
% 
% plt = figure
% s = scatter(ensembleCost(5,:), ensembleCost(6,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:1475);
%  
% plt = figure
% s = scatter(ensembleCost(1,:), ensembleCost(3,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:1475);
% 
% 
% % ich brauche aber eigene cost functions: likelihood für scalar und IQD für
% % vector data
% 
% % -> werden in costFunction2 eingebracht
% % damit costFunction2 für FullSpectrum funkctioniert, muss modData auch
% % .sizeFull haben -> in matchModOutput2Data2.m angepasst
% 
% 
% 
% 
% 
% 
% % interessanter wäre es, zwei metriken gegeneinander aufzutragen die nicht
% % korrelieren. also biogeochemie auf die y, ökologie auf die x achse. 
% % dafür habe ich die costfunction manipuliert. 
% plt = figure
% s = scatter(ensembleCost(8,:), ensembleCost(9,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
% xlabel(['Metric for biogeochemical tracers (' costFunctionType8 ')' ]);
% ylabel(['Metric for ecological state (' costFunctionType9 ')']); 
% title('Model fit to data for 2000 ensemble members')
% 
% hold on

% get cost of default/reference (optimal?) params

    % run model on initial (default) and fitted (reference) params 
    % get reference params as well
    parsFitted = matfile('results/multiplePredatorClasses/fittedParameters_RMS_Hellinger_ZPratio_Atlantic_final', 'Writable', true);
    parsFitted = parsFitted.Params;
    
    
    
    v0_initial = initialiseVariables(FixedParams, Params, Forc);
    v0_fitted = initialiseVariables(FixedParams, parsFitted, Forc);
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
    [out_initial, auxVars_initial] = integrateTrajectories(FixedParams, Params, Forc, v0_initial, ... 
        FixedParams.odeIntegrator, FixedParams.odeOptions);
    toc
    
    tic; disp('.. started at'); disp(datetime('now'))
    [out_fitted, auxVars_fitted] = integrateTrajectories(FixedParams, parsFitted, Forc, v0_fitted, ... 
        FixedParams.odeIntegrator, FixedParams.odeOptions);
    toc

    % fit model output to obs data
%     modData = matchModOutput2Data(out, auxVars, Data, FixedParams, ...
%         'fitToFullSizeSpectra', false); % only .scalar and .size
%     modData2 = matchModOutput2Data(out, auxVars, Data, FixedParams, ...
%         'fitToFullSizeSpectra', true); % output is only .scalar; no .size or .sizeFull -> true does not work. 
%     
    modData3_initial = matchModOutput2Data2(out_initial, auxVars_initial, Data, FixedParams, ...
        'fitToFullSizeSpectra', true);
    modData3_fitted = matchModOutput2Data2(out_fitted, auxVars_fitted, Data, FixedParams, ...
        'fitToFullSizeSpectra', true);
    % cellDensity and BiovolumeDensity are now extracted for full size 
    % spectra of 10 trajectories per sampling event. for costfunctions 
    % manipulated from aidans code to work, they cost values need to be 
    % averaged first. 
    

    
    % calculate cost
%     [referenceCost8, referenceCostComponents8] = costFunction2('label', costFunctionType8, ...
%             'Data', Data, 'modData', modData);
%     [referenceCost9, referenceCostComponents9] = costFunction2('label', costFunctionType9, ...
%             'Data', Data, 'modData', modData);
%         
        % also IQD and likelihood
    [referenceCostIQD_initial, referenceCostComponentsIQD_initial] = costFunction2('label', 'IQD_vectorFullSpectrum', ...
    'Data', Data, 'modData', modData3_initial); 
    [referenceCostL_initial, referenceCostComponentsL_initial] = costFunction2('label', 'syntheticLikelihood_scalar', ...
    'Data', Data, 'modData', modData3_initial); 
    [referenceCostHD_initial, referenceCostComponentsHD_initial] = costFunction2('label', 'Hellinger_vectorFullSpectrum', ...
    'Data', Data, 'modData', modData3_initial); 

    [referenceCostIQD_fitted, referenceCostComponentsIQD_fitted] = costFunction2('label', 'IQD_vectorFullSpectrum', ...
    'Data', Data, 'modData', modData3_fitted); 
    [referenceCostL_fitted, referenceCostComponentsL_fitted] = costFunction2('label', 'syntheticLikelihood_scalar', ...
    'Data', Data, 'modData', modData3_fitted); 
    [referenceCostHD_fitted, referenceCostComponentsHD_fitted] = costFunction2('label', 'Hellinger_vectorFullSpectrum', ...
    'Data', Data, 'modData', modData3_fitted); 
        
        % test if I can get fullSpectrum cost with aidans functions
        %modData3.size = modData3.sizeFull; 
        
%         [costFS costComponentsFS] = costFunction2('label', 'RMS_HellingerFullSpectrum' ,... 
%             'Data', Data, 'modData', modData3);
%         % works, with some modifications
        % cost of 'Hellinger_groupWaterOrigin_vector' and
        % 'RMS_HellingerFullSpectrum' are similar but not equal! (this is
        % good)
        

        
%     
% scatter(referenceCost8, referenceCost9,  'r')
% scatter(referenceCost8, costFS,  'r')
% scatter(referenceCostIQD, referenceCostL,  'r')
% hold off


% summer
plt = figure
s = scatter(ensembleCost.summer.IQD, ensembleCost.summer.L)
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (IQD\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for 2000 ensemble members, summer config')

% hold on
% scatter(referenceCostIQD_initial, referenceCostL_initial,  'r')
% scatter(referenceCostIQD_fitted, referenceCostL_fitted,  'g')
% hold off

% autumn
plt = figure
s = scatter(ensembleCost.autumn.IQD, ensembleCost.autumn.L)
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (IQD\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for 2000 ensemble members, autumn config')


% den gleichen plot wiederholen, aber auf der x-Achse hellinger distance
% (FullSpectrum) statt IQD
plt = figure
s = scatter(ensembleCost.summer.HD, ensembleCost.summer.L)
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (Hellinger\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for 2000 ensemble members')

hold on
scatter(referenceCostHD_initial, referenceCostL_initial,  'r')
scatter(referenceCostHD_fitted, referenceCostL_fitted,  'g')
hold off

% summer und autumn plots mit farbe fur ensMember

% summer
plt = figure
s = scatter(ensembleCost.summer.IQD, ensembleCost.summer.L, [], linspace(1,10,2000))
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (IQD\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for 2000 ensemble members, summer config')
% autumn
plt = figure
s = scatter(ensembleCost.autumn.IQD, ensembleCost.autumn.L, [], linspace(1,10,2000))
s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
xlabel(['Metric for ecological state (IQD\_vectorFullSpectrum)' ]);
ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
title('Model fit to data for 2000 ensemble members, autumn config')


% 
% 
% 
% 
% % get ensemble cost from old functions; type 9 is 'Hellinger_groupWaterOrigin_vector'
% ensembleCostOld = matfile([Directories.resultsDir 'ensembleCost_binned.m']);
% ensembleCostOld = ensembleCostOld.cost;
% ensembleCostHellinger = ensembleCostOld(9,:);
% 
% % HD (grouped) für initial and fitted 
% 
% modData_initial = matchModOutput2Data(out_initial, auxVars_initial, Data, FixedParams, ...
%     'fitToFullSizeSpectra', false);
% modData_fitted = matchModOutput2Data(out_fitted, auxVars_fitted, Data, FixedParams, ...
%     'fitToFullSizeSpectra', false);
% 
% [referenceCostHD_initial, referenceCostComponentsHD_initial] = costFunction2('label', 'Hellinger_groupWaterOrigin_vector', ...
%         'Data', Data, 'modData', modData_initial);
% [referenceCostHD_fitted, referenceCostComponentsHD_fitted] = costFunction2('label', 'Hellinger_groupWaterOrigin_vector', ...
%     'Data', Data, 'modData', modData_fitted);
% 
% plt = figure
% s = scatter(ensembleCostHellinger, ensembleCost(2,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
% xlabel(['Metric for ecological state (Hellinger\_groupWaterOrigin\_vector)' ]);
% ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
% title('Model fit to data for 2000 ensemble members')
% 
% hold on
% scatter(referenceCostHD_initial, referenceCostL_initial,  'r')
% scatter(referenceCostHD_fitted, referenceCostL_fitted,  'g')
% hold off
% 
% 
% % und HD UNGROUPED
% % get ensemble cost from old functions; type 10 is 'Hellinger_vector'
% ensembleCostHellingerUG = ensembleCostOld(10,:);
% 
% % HD (ungrouped) für initial and fitted 
% [referenceCostHD_UG_initial, referenceCostComponentsHD_UG_initial] = costFunction2('label', 'Hellinger_vector', ...
%         'Data', Data, 'modData', modData_initial);
% [referenceCostHD_UG_fitted, referenceCostComponentsHD_UG_fitted] = costFunction2('label', 'Hellinger_vector', ...
%     'Data', Data, 'modData', modData_fitted);
% 
% plt = figure
% s = scatter(ensembleCostHellingerUG, ensembleCost(2,:))
% s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Member ID', 1:2000);
% xlabel(['Metric for ecological state (Hellinger\_vector)' ]);
% ylabel(['Metric for biogeochemical tracers (syntheticLikelihood\_scalar)']); 
% title('Model fit to data for 2000 ensemble members')
% 
% hold on
% scatter(referenceCostHD_UG_initial, referenceCostL_initial,  'r')
% scatter(referenceCostHD_UG_fitted, referenceCostL_fitted,  'g')
% hold off



