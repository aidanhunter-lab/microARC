%% Size-structured 1D NPZD model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Generate model outputs from suite of parameter sets
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Refresh the workspace
clear; close all; delete(gcp('nocreate')); clc
% Include all subdirectories within search path
addpath(genpath(fileparts(which('fit_parameters'))))

rng(1) % set random seed


%% Set up model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Store folders/filenames of data and saved parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Set parFile = [] for default initial parameter values (hard-coded, not 
% loaded) or parFile = 'filename.mat' to use saved values.
% fileSuffix = 'meanCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot';
fileSuffix = 'meanCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot_removeParams';
% fileSuffix = 'smoothCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot_removeParams';
parFile = ['parameterInitialValues_' fileSuffix, '.mat'];
Directories = setDirectories('bioModel', 'multiplePredatorClasses', ...
    'parFile', parFile);
display(Directories)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Load optimisation results or use default set-up
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

loadPriorOptimisation = true; % To use default set-up specify as false

switch loadPriorOptimisation
    
    case true
        resultsFile = ['fittedParameters_' fileSuffix, '.mat'];
        % Load stored results
        [~, results, ~, ~, boundsLower, boundsUpper, Data, Forc, FixedParams, Params, v0] = ...
            loadOptimisationRun(resultsFile);
        Params = updateParameters(Params, FixedParams, results.optPar);
        
        
    case false
        
        % USE THE DEFAULT SET-UP HERE... BUT WAIT UNTIL I'VE SORTED OUT THE
        % LOADED OPTION BEFORE SETTING THIS UP, SO THAT I KNOW WHAT'S NEEDED
        
        
        
%         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         % Create main structs needed to run model: Forc, FixedParams, Params, Data
%         %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%         
%         % Use only a single trajectory per sampling event, but note that samples
%         % taken near Arctic-Atlantic water boundaries are not easily descernable
%         % from forcing data when single trajectories are used.
%         numTraj = 1;
%         % Model set-up may be modified by passing some extra arguments as
%         % name-value pairs (preferable), or directly modified within modelSetUp.m.
%         % Useful name-value pairs for modelSetUp.m include: useTraj, ESDmin, ESDmax, nsizes
%         [Forc, FixedParams, Params, Data] = modelSetUp(Directories, ...
%             'displayAllOutputs', true, 'numTraj', numTraj);
%         

        
end


%% Sensitivity analysis

% Evaluate the sensitivity of the (components of the) cost function to each
% fitted parameter.

% Vary the values of all tunable parameters while holding all other
% parameters fixed. Evaluate the cost for each value then make plots to
% assess which parameters have most influence upon model fit to data, and
% which may be reasonable to omit from optimisation...

parRes = 50;
parNames = FixedParams.tunePars;
npars = length(parNames);
parVals = nan(parRes, npars);

for i = 1:npars
    parVals(:,i) = linspace(boundsLower(i), boundsUpper(i), parRes);
end

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end


% FixedParams.fitToFullSizeSpectra = false;


ParamSensitivity.parNames = FixedParams.tunePars;
ParamSensitivity.value = parVals;
% ParamSensitivity.cost = nan(size(parVals));
% ParamSensitivity.costComponentFields = [];
% ParamSensitivity.costComponents = nan(size(parVals));

% For each parameter set, store the model output required to calculate the
% cost function.

tic; disp('.. started at'); disp(datetime('now'))
progress = waitbar(0, 'progress');
for i = 1:npars
    parName = FixedParams.tunePars{i};
    x = cellfun(@(z) Params.(z), FixedParams.tunePars);
    for j = 1:parRes
        x(i) = parVals(j,i);
        % Rather than calling the 'costCalc' wrapper function, this could
        % be deconstructed to integrate the model then, separately, to call
        % multiple different cost function types...
        % DO THIS. SAVE THE MODELLED OUTPUT, THEN CALL DIFFERENT COST
        % FUNCTIONS TO PLOT SENSIVITIES FOR EACH. ONLY PROBLEM IS THAT
        % PARAMS ARE HELD FIXED AT OPTIMAL VALUES DETERMINED BY THE COST
        % FUNCTION, SO USING SAME MODELLED VALUES FOR EACH PLOT IS
        % INAPPROPRIATE, BUT THE PLOTS WILL STILL BE USEFUL...
        [modData, out, auxVars] = ensembleMakeOutput(...
            x, FixedParams, Params, Forc, Data, v0, ...
            FixedParams.odeIntegrator, FixedParams.odeOptions, ...
            'returnExtra', 'all');
        if j > 1 || i > 1
            ParamSensitivity.modData.scalar.Value(:,j,i) = modData.scalar.Value;
            if isfield(modData.scalar, 'scaled_Value')
                ParamSensitivity.modData.scalar.scaled_Value(:,j,i) = modData.scalar.scaled_Value;
            end
            ParamSensitivity.modData.size.Value(:,:,j,i) = modData.size.Value;
        elseif i == 1 && j == 1
            ParamSensitivity.modData = modData;
            % Expand dimensions to store output over parameter ensemble
            ParamSensitivity.modData.scalar.Value = nan(size(modData.scalar.Value, 1), parRes, npars);
            ParamSensitivity.modData.scalar.Value(:,j,i) = modData.scalar.Value;
            if isfield(modData.scalar, 'scaled_Value')
                ParamSensitivity.modData.scalar.scaled_Value = nan(size(modData.scalar.scaled_Value, 1), parRes, npars);
                ParamSensitivity.modData.scalar.scaled_Value(:,j,i) = modData.scalar.scaled_Value;
            end
            ParamSensitivity.modData.size.Value = nan(size(modData.size.Value, 1), size(modData.size.Value, 2), parRes, npars);
            ParamSensitivity.modData.size.Value(:,:,j,i) = modData.size.Value;
        end
%         [cost, costComponents, modData, out, auxVars] = costCalc(...
%             x, FixedParams, Params, Forc, Data, v0, ...
%             FixedParams.odeIntegrator, FixedParams.odeOptions, ... 
%             'selectFunction', FixedParams.costFunction, ...
%             'returnExtra', 'all');
%         if i == 1 && j == 1
%             costComponentFields = fieldnames(costComponents);
%             nFields = length(costComponentFields);
%             ParamSensitivity.costComponentFields = costComponentFields;
%             ParamSensitivity.costComponents = repmat(reshape( ... 
%                 ParamSensitivity.costComponents, ... 
%                 [1, size(ParamSensitivity.costComponents)]), [nFields, 1, 1]);
%             clear costComponentFields nFields
%         end
%         ParamSensitivity.cost(j,i) = cost;
%         ParamSensitivity.costComponents(:,j,i) = struct2array(costComponents);
        waitbar(((i-1) * parRes + j) / (npars*parRes), progress)
    end
end
Time = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
disp(['Run-time: ' num2str(floor(Time)) ' hrs, ' ...
    num2str(floor(mod(60*Time,60))) ' mins'])



% Sort out how to save this stuff...

fileName = ['paramSensitivity_' fileSuffix];

m = matfile(fileName, 'Writable', true);
% m.results = ParamSensitivity;
ParamSensitivity = m.results;


% switch saveParams, case true
%     % Fitted parameters are saved as the 'results' struct, and the 
%     % associated model set-up is saved within each subsequent struct
%     if ~exist('v0', 'var'), v0 = []; end
%     saveOptimisationRun(fileName_results, results, Data, Forc, FixedParams, Params, v0);
% end
% 


% Calculate cost function values given various parameter sets using model
% output stored in ParamSensitivity

% ParamSensitivity.cost = nan(size(parVals));
% ParamSensitivity.costComponentFields = [];
% ParamSensitivity.costComponents = nan(size(parVals));
% 

selectFunction = 'LeastAbsErr_Hellinger';

for i = 1:npars
    parName = ParamSensitivity.parNames{i};
    x = cellfun(@(z) Params.(z), FixedParams.tunePars);
    for j = 1:parRes
        x(i) = ParamSensitivity.value(j,i);
        
        modData = ParamSensitivity.modData;
        modData.scalar.Value = modData.scalar.Value(:,j,i);
        if isfield(modData.scalar,'scaled_Value')
            modData.scalar.scaled_Value = modData.scalar.scaled_Value(:,j,i);
        end
        modData.size.Value = modData.size.Value(:,:,j,i);
        [cost, costComponents] = costFunction('label', selectFunction, ...
            'Data', Data, 'modData', modData);
        if i > 1 && j > 1
            ParamSensitivity.cost(j,i) = cost;
            ParamSensitivity.costComponents(j,i,:) = struct2array(costComponents);
        elseif i == 1 && j == 1
            ParamSensitivity.cost = nan(size(ParamSensitivity.value));
            ParamSensitivity.cost(j,i) = cost;
            costFields = fieldnames(costComponents);
            ParamSensitivity.costFields = costFields;
            ParamSensitivity.costComponents = nan([size(ParamSensitivity.cost), length(costFields)]);
            ParamSensitivity.costComponents(j,i,:) = struct2array(costComponents);
        end
    end
end





% Plots of conditional sensitivities (each parameter's sensitivity
% conditioned on fixed/fitted values of all other parameters)

parNames
negExpTransform = false(npars,1);
negExpTransform(ismember(parNames, {'Gmax_b','pmax_b','Qmax_delQ_b','aN_QC_b'})) = true;

plt1 = figure;
plt1.Units = 'inches';
plt1.Position = [0 0 16 8];

nrows = floor(sqrt(npars));
ncols = ceil(npars / nrows);

for i = 1:npars
    subplot(nrows, ncols, i)
    p = parNames{i};
    ind = strcmp(ParamSensitivity.parNames, p);
    negExp = negExpTransform(i);
    if ~negExp
        value = ParamSensitivity.value(:,ind);
        valueFix = Params.(p);
    else
        value = -exp(ParamSensitivity.value(:,ind));
        valueFix = -exp(Params.(p));
    end
    cost = ParamSensitivity.cost(:,ind);
    scatter(value, cost)
    hold on
    scatter(valueFix, interp1(value, cost, valueFix), ...
        'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
    hold off
    xlabel(p, 'Interpreter', 'none')
    ylabel('cost')
end


plt2 = figure;
plt2.Units = 'inches';
plt2.Position = [0 0 16 8];

nrows = floor(sqrt(npars));
ncols = ceil(npars / nrows);

costOmit = contains(ParamSensitivity.costComponentFields, 'Tot');
costComponents = ParamSensitivity.costComponents(~costOmit,:,:);

for i = 1:npars
    subplot(nrows, ncols, i)
    p = parNames{i};
    ind = strcmp(ParamSensitivity.parNames, p);
    negExp = negExpTransform(i);
    if ~negExp
        value = ParamSensitivity.value(:,ind);
        valueFix = Params.(p);
    else
        value = -exp(ParamSensitivity.value(:,ind));
        valueFix = -exp(Params.(p));
    end
    cost = costComponents(:,:,ind);
    plot(value, cost)
    gc = gca;
    hold on
    scatter(valueFix, gc.YLim(1), 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
    hold off
    xlabel(p, 'Interpreter', 'none')
    ylabel('cost')
end
i = i + 2;
sp = subplot(nrows, ncols, i); plot(1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan); set(sp, 'Visible', 'off');

legLabel = ParamSensitivity.costComponentFields(~costOmit);


legend(sp, legLabel,'Location', 'west', 'Interpreter', 'none');



plt3 = figure;
plt3.Units = 'inches';
plt3.Position = [0 0 16 8];

nrows = floor(sqrt(npars));
ncols = ceil(npars / nrows);

costOmit = contains(ParamSensitivity.costComponentFields, 'Tot');
costComponents = ParamSensitivity.costComponents(~costOmit,:,:);

scalarInd = contains(ParamSensitivity.costComponentFields(~costOmit), Data.scalar.obsInCostFunction);
sizeInd = contains(ParamSensitivity.costComponentFields(~costOmit), Data.size.obsInCostFunction);

for i = 1:npars
    subplot(nrows, ncols, i)
    p = parNames{i};
    ind = strcmp(ParamSensitivity.parNames, p);
    negExp = negExpTransform(i);
    if ~negExp
        value = ParamSensitivity.value(:,ind);
        valueFix = Params.(p);
    else
        value = -exp(ParamSensitivity.value(:,ind));
        valueFix = -exp(Params.(p));
    end
    cost = costComponents(:,:,ind);
    
    costScalar = mean(cost(scalarInd,:));
    costSize = mean(cost(sizeInd,:));
   
    plot(value, costScalar)
    gc = gca;
    hold on
    plot(value, costSize)
    
    scatter(valueFix, gc.YLim(1), 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
    hold off
    xlabel(p, 'Interpreter', 'none')
    ylabel('cost')
end
i = i + 1;
sp = subplot(nrows, ncols, i); plot(1, nan, 1, nan); set(sp, 'Visible', 'off');

legLabel = {'scalar','size'};

legend(sp, legLabel,'Location', 'west', 'Interpreter', 'none');



% Plots of marginal sensitivities -- cannot plot the cost components
% without rerunning the model and explicitly storing those costs because
% the ga output function can only account for a scalar total cost...

popHist = results.populationHistory;
costHist = results.scoreHistory;


plt4 = figure;
plt4.Units = 'inches';
plt4.Position = [0 0 16 8];

nrows = floor(sqrt(npars));
ncols = ceil(npars / nrows);

for i = 1:npars
    subplot(nrows, ncols, i)
    p = parNames{i};
    ind = strcmp(results.parNames, p);
    negExp = negExpTransform(i);
    parVals = squeeze(popHist(:,ind,:));
    if negExp
        parVals = -exp(parVals);
    end
    scatter(parVals(:), costHist(:))
    xlabel(p, 'Interpreter', 'none')
    ylabel('cost')
end



%% Try different cost functions for the scalar (nutrient) data

% Generate model output for an ensemble of parameter sets, varying only
% those few parameters that are not informed by the size data and entirely
% informed by the scalar data.

parNames = FixedParams.tunePars;
vPars = {'wPOM1', 'rDON', 'rPON'};

parRes = 15;
npars = length(vPars);
nVals = parRes ^ npars;
parVals = nan(parRes, npars);

for i = 1:npars
    j = strcmp(parNames, vPars{i});
    parVals(:,i) = linspace(boundsLower(j), boundsUpper(j), parRes);
end


ParamEnsemble.parNames = vPars;
ParamEnsemble.value = parVals;
[X1, X2, X3] = ndgrid(parVals(:,1), parVals(:,2), parVals(:,3));
ParamEnsemble.allValues = cat(2, X1(:), X2(:), X3(:));

tic; disp('.. started at'); disp(datetime('now'))
progress = waitbar(0, 'progress');
x = cellfun(@(z) Params.(z), FixedParams.tunePars);
for i = 1:nVals
    xi = ParamEnsemble.allValues(i,:);
    for j = 1:npars
        x(strcmp(parNames, vPars{j})) = xi(j);
    end
%     x(contains(parNames, vPars)) = xi;
    [modData, out, auxVars] = ensembleMakeOutput(...
        x, FixedParams, Params, Forc, Data, v0, ...
        FixedParams.odeIntegrator, FixedParams.odeOptions, ...
        'returnExtra', 'all');
    if i > 1
        ParamEnsemble.modData.scalar.Value(:,i) = modData.scalar.Value;
        if isfield(modData.scalar, 'scaled_Value')
            ParamEnsemble.modData.scalar.scaled_Value(:,i) = modData.scalar.scaled_Value;
        end
        ParamEnsemble.modData.size.Value(:,:,i) = modData.size.Value;
    else
        ParamEnsemble.modData = modData;
        % Expand dimensions to store output over parameter ensemble
        ParamEnsemble.modData.scalar.Value = nan(size(modData.scalar.Value, 1), nVals);
        ParamEnsemble.modData.scalar.Value(:,1) = modData.scalar.Value;
        if isfield(modData.scalar, 'scaled_Value')
            ParamEnsemble.modData.scalar.scaled_Value = nan(size(modData.scalar.scaled_Value, 1), nVals);
            ParamEnsemble.modData.scalar.scaled_Value(:,1) = modData.scalar.scaled_Value;
        end
        ParamEnsemble.modData.size.Value = nan(size(modData.size.Value, 1), size(modData.size.Value, 2), nVals);
        ParamEnsemble.modData.size.Value(:,:,1) = modData.size.Value;
    end
    waitbar(i/nVals, progress)
end
Time = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
disp(['Run-time: ' num2str(floor(Time)) ' hrs, ' ...
    num2str(floor(mod(60*Time,60))) ' mins'])
    

% Save
fileSuffix = 'meanCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot_removeParams';
fileName = ['paramEnsemble' fileSuffix];

m2 = matfile(fileName, 'Writable', true);
% m2.results = ParamEnsemble;
ParamEnsemble = m2.results;


% Use ensemble model output to assess nutrient component of cost -- using
% various different functions.

% modData_ = ParamEnsemble.modData;
% 
% modData.scalar = modData_.scalar;
% modData.scalar.Value = modData.scalar.Value(:,1);
% modData.scalar.scaled_Value = modData.scalar.scaled_Value(:,1);
% modData.size = modData_.size;
% modData.size.Value = modData.size.Value(:,:,1);

selectFunction = 'meanCDFdist_Hellinger';
selectFunction = 'smoothCDFdist_Hellinger';
selectFunction = 'LeastAbsErr_Hellinger';

% [cost1, costComponents1] = costFunction('label', selectFunction, ...
%     'Data', Data, 'modData', modData);
% [cost2, costComponents2] = costFunction('label', selectFunction, ...
%     'Data', Data, 'modData', modData);


progress = waitbar(0, 'progress');
for i = 1:nVals
    modData.scalar = ParamEnsemble.modData.scalar;
    modData.size = ParamEnsemble.modData.size;
    modData.scalar.Value = modData.scalar.Value(:,i);
    modData.scalar.scaled_Value = modData.scalar.scaled_Value(:,i);
    modData.size.Value = modData.size.Value(:,:,i);
    [cost, costComponents] = costFunction('label', selectFunction, ...
        'Data', Data, 'modData', modData);
    if i > 1
        ParamEnsemble.cost(i) = cost;
        ParamEnsemble.costComponents(i,:) = struct2array(costComponents);
    else
        ParamEnsemble.cost = nan(size(ParamEnsemble.allValues, 1), 1);
        ParamEnsemble.cost(i) = cost;
        costFields = fieldnames(costComponents);
        ParamEnsemble.costFields = costFields;
        ParamEnsemble.costComponents = nan(size(ParamEnsemble.cost, 1), length(costFields));
        ParamEnsemble.costComponents(i,:) = struct2array(costComponents);
    end
    waitbar(i/nVals, progress)
end


% Make surface or contour plots to display cost function values relative to
% two parameters -- slices out of the overall cost surface.

allValues = ParamEnsemble.allValues;
allValues = reshape(allValues, [parRes, parRes, parRes, 3]);
cost = ParamEnsemble.cost;
cost = reshape(cost, [parRes, parRes, parRes]);
costComponents = ParamEnsemble.costComponents;
costComponents = reshape(costComponents, [parRes, parRes, parRes, length(ParamEnsemble.costFields)]);

costTot = cost;
costNut = mean(costComponents(:,:,:,1:4), 4);
costSize = mean(costComponents(:,:,:,[5,7]), 4);

plt = figure;
plt.Units = 'inches';
plt.Position = [0 0 16 16];

Zmin = nan(parRes, 3);
for i = 1:parRes
    ziTot = costTot(i,:,:);
    ziNut = costNut(i,:,:);
    ziSize = costSize(i,:,:);
    Zmin(i,1) = min(ziNut(:));
    Zmin(i,2) = min(ziSize(:));
    Zmin(i,3) = min(ziTot(:));
%     Zmin(i) = min(zi(:));
end
pi1 = Zmin(:,1) == min(Zmin(:,1)); % index best-fitting slice
pi2 = Zmin(:,2) == min(Zmin(:,2)); % index best-fitting slice
pi3 = Zmin(:,3) == min(Zmin(:,3)); % index best-fitting slice

pi=8;pi1=pi;pi2=pi;pi3=pi;

Ztot = squeeze(costTot(pi3,:,:));
Znut = squeeze(costNut(pi1,:,:));
Zsize = squeeze(costSize(pi2,:,:));
x = ParamEnsemble.value(:,2);
y = ParamEnsemble.value(:,3);

subplot(3,3,1)
surf(x,y,Znut)
xlabel(vPars{2})
ylabel(vPars{3})
title('nutrient')

subplot(3,3,2)
surf(x,y,Zsize)
xlabel(vPars{2})
ylabel(vPars{3})
title('size')

subplot(3,3,3)
surf(x,y,Ztot)
xlabel(vPars{2})
ylabel(vPars{3})
title('total')


Zmin = nan(parRes, 3);
for i = 1:parRes
    ziTot = costTot(:,i,:);
    ziNut = costNut(:,i,:);
    ziSize = costSize(:,i,:);
    Zmin(i,1) = min(ziNut(:));
    Zmin(i,2) = min(ziSize(:));
    Zmin(i,3) = min(ziTot(:));
%     Zmin(i) = min(zi(:));
end
pi1 = Zmin(:,1) == min(Zmin(:,1)); % index best-fitting slice
pi2 = Zmin(:,2) == min(Zmin(:,2)); % index best-fitting slice
pi3 = Zmin(:,3) == min(Zmin(:,3)); % index best-fitting slice

Ztot = squeeze(costTot(:,pi3,:));
Znut = squeeze(costNut(:,pi1,:));
Zsize = squeeze(costSize(:,pi2,:));
x = ParamEnsemble.value(:,1);
y = ParamEnsemble.value(:,3);

subplot(3,3,4)
surf(x,y,Znut)
xlabel(vPars{1})
ylabel(vPars{3})
title('nutrient')

subplot(3,3,5)
surf(x,y,Zsize)
xlabel(vPars{1})
ylabel(vPars{3})
title('size')

subplot(3,3,6)
surf(x,y,Ztot)
xlabel(vPars{1})
ylabel(vPars{3})
title('total')

Zmin = nan(parRes, 3);
for i = 1:parRes
    ziTot = costTot(:,:,i);
    ziNut = costNut(:,:,i);
    ziSize = costSize(:,:,i);
    Zmin(i,1) = min(ziNut(:));
    Zmin(i,2) = min(ziSize(:));
    Zmin(i,3) = min(ziTot(:));
%     Zmin(i) = min(zi(:));
end
pi1 = Zmin(:,1) == min(Zmin(:,1)); % index best-fitting slice
pi2 = Zmin(:,2) == min(Zmin(:,2)); % index best-fitting slice
pi3 = Zmin(:,3) == min(Zmin(:,3)); % index best-fitting slice

Ztot = squeeze(costTot(:,:,pi3));
Znut = squeeze(costNut(:,:,pi1));
Zsize = squeeze(costSize(:,:,pi2));
x = ParamEnsemble.value(:,1);
y = ParamEnsemble.value(:,2);

subplot(3,3,7)
surf(x,y,Znut)
xlabel(vPars{1})
ylabel(vPars{2})
title('nutrient')

subplot(3,3,8)
surf(x,y,Zsize)
xlabel(vPars{1})
ylabel(vPars{2})
title('size')

subplot(3,3,9)
surf(x,y,Ztot)
xlabel(vPars{1})
ylabel(vPars{2})
title('total')









Zmin = nan(parRes, 1);
for i = 1:parRes
    zi = costComponents(:,:,i);
    Zmin(i) = min(zi(:));
end
pi = Zmin == min(Zmin); % index best-fitting slice

Z = squeeze(costComponents(:,:,pi));
x = ParamEnsemble.value(:,1);
y = ParamEnsemble.value(:,2);
surf(x,y,Z)
xlabel(vPars{1})
ylabel(vPars{2})


Zmin = nan(parRes, 1);
for i = 1:parRes
    zi = costComponents(i,:,:);
    Zmin(i) = min(zi(:));
end
pi = Zmin == min(Zmin); % index best-fitting slice

Z = squeeze(costComponents(8,:,:));
x = ParamEnsemble.value(:,2);
y = ParamEnsemble.value(:,3);
surf(x,y,Z)
xlabel(vPars{2})
ylabel(vPars{3})






Z = squeeze(costComponents(:,1,:));

x = ParamEnsemble.value(:,2);
y = ParamEnsemble.value(:,3);

surf(x,y,Z)
xlabel(vPars{2})
ylabel(vPars{3})


size(costComponents)
size(allValues)



figure
x = ParamEnsemble.allValues(:,1);
y = mean(ParamEnsemble.costComponents(:,1:4), 2);
scatter(x,y)

figure
x = ParamEnsemble.allValues(:,2);
y = mean(ParamEnsemble.costComponents(:,1:4), 2);
scatter(x,y)

figure
x = ParamEnsemble.allValues(:,3);
y = mean(ParamEnsemble.costComponents(:,1:4), 2);
scatter(x,y)



figure

allValues = ParamEnsemble.allValues;
allValues = reshape(allValues, [parRes, parRes, parRes, 3]);
cost = ParamEnsemble.cost;
cost = reshape(cost, [parRes, parRes, parRes]);
costComponents = ParamEnsemble.costComponents;
costComponents = reshape(costComponents, [parRes, parRes, parRes, length(ParamEnsemble.costFields)]);

p = ParamEnsemble.allValues(:,1);
c = ParamEnsemble.cost;

j = 1;
p = allValues(j,:,:,3);
% c = cost(:,:,j);
c1 = mean(costComponents(j,:,:,1:4), 4); % nutrient data
c2 = mean(costComponents(j,:,:,[5, 7]), 4); % size data


scatter(p(:),c1(:), 'MarkerEdgeColor', [0 0 1])
hold on
scatter(p(:),c2(:), 'MarkerEdgeColor', [1 0 0])

for j = 2:parRes
    p = allValues(:,:,j,3);
    c1 = mean(costComponents(:,:,j,1:4), 4); % nutrient data
    c2 = mean(costComponents(:,:,j,[5, 7]), 4); % size data
    scatter(p(:),c1(:), 'MarkerEdgeColor', [0 0 1])
    scatter(p(:),c2(:), 'MarkerEdgeColor', [1 0 0])
end

hold off



hold on
for j = 2:parRes
    p = allValues(:,:,j,1);
    c = cost(:,:,j);
    scatter(p(:),c(:))
end

hold off



hist(cost(:))

xx = costComponents(:,:,:,1:4);
xx = costComponents(:,:,:,[5,7]);
xx = mean(xx, 4);
hist(xx(:))



%%

% Plot params against each other to identify possible correlations
ind1 = strcmp(parNames, 'wPOM1');
ind2 = strcmp(parNames, 'rDON');
ind3 = strcmp(parNames, 'rPON');

ind1 = strcmp(parNames, 'Vmax_QC_a');
ind2 = strcmp(parNames, 'Vmax_QC_b');
ind3 = strcmp(parNames, 'aN_QC_a');
ind4 = strcmp(parNames, 'aN_QC_b');

parVals1 = squeeze(popHist(:,ind1,1:end));
parVals2 = squeeze(popHist(:,ind2,1:end));
parVals3 = squeeze(popHist(:,ind3,1:end));
parVals4 = squeeze(popHist(:,ind4,1:end));

figure
scatter(parVals1(:), parVals2(:))
scatter(parVals1(:), parVals3(:))
scatter(parVals2(:), parVals3(:))


figure
histogram(costHist(:,end))

plot(min(costHist))


%% Generate parameters

% Choose which parameters to vary and the cost function(s) to evaluate.

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NOTE...
% To compute in reasonable time, the number of parameters selected to vary
% must be low -- otherwise the "curse of dimensionality" strikes!
% We should therefore select small groups of parameters to vary in batches
% while holding other parameters at constant values. If we choose, say, 4
% parameters to vary and specify 8 separate values for each one, then 
% generating all model outputs will require 8^4=4096 function evaluations.
% As each function evaluation is computationally costly, using larger
% parameter batches and/or higher parameter resolutions will quickly become
% unfeasible...
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

costFunctionType = 'meanCDFdist_Hellinger';
fitTrajectories = []; % May be set as 'Atlantic' or 'Arctic', but probably best left as [] when generating ensemble outputs.
fitToFullSizeSpectra = false; % Set to true if using cost function that compares model output to full (unbinned) size spectra data.
parSuite = {'Vmax_QC_b','aN_QC_b','pmax_b','Gmax_b'}; % Parameters to vary
parRes = 8;

[FixedParams, Params, Forc, Data] = ... 
    ensembleOptions(FixedParams, Params, Forc, Data, ...
    'costFunctionType', costFunctionType, ...
    'fitTrajectories', fitTrajectories, ...
    'fitToFullSizeSpectra', fitToFullSizeSpectra, ...
    'parSuite', parSuite, ...
    'parRes', parRes);

% Create a suite of parameters sets over which to evaluate the model.
% Use a latin hypercube to generate evenly distributed parameters.
np = 8; % Number of values per parameter

nf = np .^ 4;

lhsdesign()



% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
v0 = initialiseVariables(FixedParams, Params, Forc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [size, depth, nutrient]
%                      organic matter      [type, depth, nutrient]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









% Choose which parameters to tune, the cost function, and numerical tuning
% algorithm. Can also select which data to fit / trajectories to run as
% Arctic and/or Atlantic.
niter = 50;
costFunctionType = 'meanCDFdist_Hellinger';
% costFunctionType = 'meanCDFdist_HellingerFullSpectrum';
% costFunctionType = 'meanCDFdist_HellingerFullSpectrum_averagedEventsDepths';

fitTrajectories = 'Atlantic';
% Different cost functions may use binned size data (integrated within
% modelled size class intervals) or may fit to the full size spectra data.
% This choice affects how model outputs are extracted to match data.
fitToFullSizeSpectra = false;
[FixedParams, Params, Forc, Data] = ...
    optimisationOptions(FixedParams, Params, Forc, Data, ...
    'niter', niter, ...
    'costFunctionType', costFunctionType, ...
    'fitTrajectories', fitTrajectories, ...
    'fitToFullSizeSpectra', fitToFullSizeSpectra);

% Optional arguments (e.g. 'niter') may be included as name-value pairs,
% otherwise default values are used.
% It is important to specify 'costFunctionType' as one of the options
% available in costFunction.m.
% The 'fitTrajectories' option is also important. It is used to select sets
% of trajectories to use within the parameter optimisation. Set 
% fitTrajectories = 'Atlantic' or 'Arctic' to use trajectories originating
% from either region and the associated fitting data, or set
% fitTrajectories = 'all' to fit to all data using using trajectories
% originating from Arctic or Atlantic. Omit fitTrajectories or set it empty
% ([]) to optimise using all trajectories and fitting to size data that is
% aggregated over all samples.

% Store state variables in array v0: 1st dimension = variable
%                                    2nd dimension = location (trajectory)
v0 = initialiseVariables(FixedParams, Params, Forc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Order of variables = inorganic nutrients [depth]
%                      phytoplankton       [size, depth, nutrient]
%                      zooplankton         [size, depth, nutrient]
%                      organic matter      [type, depth, nutrient]
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

restartRun = true; % restart algorithm from a saved prior run?
switch restartRun, case true
    fileName_results = 'fittedParameters';  % saved parameters file name
%     tag = '1';                              % and identifying tag
    tag = FixedParams.costFunction;
%     tag = [tag '_Atlantic_quadraticMortality_singleTraj_relativeSizeDataOnly'];
    tag = [tag, '_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot'];
    fileName_results = fullfile(Directories.resultsDir, ...
        [fileName_results '_' tag]);
    % Load stored results    
    [~, results, ~, ~, boundsLower, boundsUpper, Data, Forc, FixedParams, Params, v0] = ...
        loadOptimisationRun(fileName_results);
    populationOld = results.populationHistory(:,:,end);
    scoresOld = results.scoreHistory(:,end);
    optimiserOptions = results.optimiserOptions;
    optimiserOptions.MaxGenerations = niter;
    optimiserOptions.InitialPopulationMatrix = populationOld;
    optimiserOptions.InitialScoresMatrix = scoresOld;
end

% Parallelise integrations over trajectories
poolObj = gcp('nocreate');
if isempty(poolObj), poolObj = parpool('SpmdEnabled', false); end

% Call optimiser
tic; disp('.. started at'); disp(datetime('now'))
[optPar, fval, exitflag, output, population, scores] = optimise(@(x) ... 
    costCalc(x, FixedParams, Params, Forc, Data, v0, FixedParams.odeIntegrator, ...
    FixedParams.odeOptions, 'selectFunction', costFunctionLabel), ... 
    npars, [], [], [], [], boundsLower, boundsUpper, [], optimiserOptions);
optimisationTime = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
disp(['Optimisation time: ' num2str(floor(optimisationTime)) ' hrs, ' ...
    num2str(floor(mod(60*optimisationTime,60))) ' mins'])

