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
% fileSuffix = 'smoothCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot_removeParams';

% fileSuffix = 'meanCDFdist_Hellinger_Atlantic_quadraticMortality_singleTraj_omitSizeDataTot_removeParams';

fileSuffix = 'RMS_Hellinger2_Atlantic_singleTraj_removeParams';

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
        
        % Make this work without loading fitted values -- default set-up...
        
end


%% Sensitivity analysis

% Evaluate the sensitivity of the (components of the) cost function to each
% fitted parameter.
% Vary the value of each tunable parameter while holding all other
% parameters fixed at some fitted (optimal) values. Store modelled output
% for each parameter set to assess parameter (conditional) sensitivities.

loadPriorRun = true;

switch loadPriorRun
    
    case true
        % Load stored results
        fileName = ['paramSensitivity_' fileSuffix];
        m = matfile(fullfile(Directories.resultsDir, fileName), 'Writable', true);
        ParamSensitivity = m.results; % load
        parNames = ParamSensitivity.parNames;
        npars = length(parNames);
        parRes = size(ParamSensitivity.value, 1);
        
    case false
        % Create ParamSensitivity struct
        parNames = FixedParams.tunePars; % parameters of interest
        npars = length(parNames); % number of parameters
        parRes = 50; % Number of values per parameter (evenly spaced between parameter bounds)
        parVals = nan(parRes, npars); % grid of values
        for i = 1:npars
            parVals(:,i) = linspace(boundsLower(i), boundsUpper(i), parRes);
        end
        % Outputs are stored in ParamSensitivity struct
        ParamSensitivity.parNames = FixedParams.tunePars;
        ParamSensitivity.value = parVals;

        % For each parameter set, store the model output, modData, required
        % to calculate the cost function. (Any model output may be stored, 
        % but due to memory restrictions we need to be selective.)
        
        tic; disp('.. started at'); disp(datetime('now'))
        progress = waitbar(0, 'progress');
        for i = 1:npars
            parName = FixedParams.tunePars{i};
            x = cellfun(@(z) Params.(z), FixedParams.tunePars);
            for j = 1:parRes
                x(i) = parVals(j,i);
                % Call wrapper function "ensembleMakeOutput" to integrate model for
                % chosen parameters and to generate modData
                [modData, out, auxVars] = ensembleMakeOutput(...
                    x, FixedParams, Params, Forc, Data, v0, ...
                    FixedParams.odeIntegrator, FixedParams.odeOptions, ...
                    'returnExtra', 'all');
                % Store model output in ParamSensitivity struct
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
                waitbar(((i-1) * parRes + j) / (npars*parRes), progress)
            end
        end
        Time = toc / 60 / 60; disp('.. finished at'); disp(datetime('now')); fprintf('\n')
        disp(['Run-time: ' num2str(floor(Time)) ' hrs, ' ...
            num2str(floor(mod(60*Time,60))) ' mins'])
        close(progress)
end

saveOutput = true;

switch saveOutput, case true
    fileName = ['paramSensitivity_' fileSuffix];
    m = matfile(fullfile(Directories.resultsDir, fileName), 'Writable', true);
    m.results = ParamSensitivity;
end


% Evaluate cost function(s) for each parameter set stored in
% ParamSensitivity to assess which parameters have most influence upon
% model fit to data, and which may be reasonable to omit from optimisation.

% Choose cost fucntion
% selectFunction = 'LeastAbsErr_Hellinger';
selectFunction = 'RMS_Hellinger2';

for i = 1:npars
    parName = ParamSensitivity.parNames{i};
    x = cellfun(@(z) Params.(z), FixedParams.tunePars); % fitted (optimal) parameters
    for j = 1:parRes
        x(i) = ParamSensitivity.value(j,i); % vary param i while params k~=i keep their fitted values
        modData = ParamSensitivity.modData;
        modData.scalar.Value = modData.scalar.Value(:,j,i); % modelled equivalents to data given selected parameters
        if isfield(modData.scalar,'scaled_Value')
            modData.scalar.scaled_Value = modData.scalar.scaled_Value(:,j,i);
        end
        modData.size.Value = modData.size.Value(:,:,j,i);
        [cost, costComponents] = costFunction('label', selectFunction, ...
            'Data', Data, 'modData', modData); % evaluate costs
        % Store in cost terms in ParamSensitivity struct
        if i > 1 || j > 1
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
negExpTransform = false(npars,1); % (I don't like the way these exponential transformed are handled... find a slicker way)
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


% Similar plot displaying cost components

plt2 = figure;
plt2.Units = 'inches';
plt2.Position = [0 0 16 8];

npars = size(ParamSensitivity.value, 2);
parRes = size(ParamSensitivity.value, 1);
costFields = ParamSensitivity.costFields;
ncost = length(costFields);

costComponents = ParamSensitivity.costComponents;
% Remove any components that evaluate to NaN
costComponents = reshape(costComponents, [npars * parRes, ncost]);
costOmit = all(isnan(costComponents));
costFields = costFields(~costOmit);
ncost = length(costFields);
costComponents = reshape(costComponents(:,~costOmit), [parRes, npars, ncost]);

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
    cost = squeeze(costComponents(:,ind,:));
    plot(value, cost)
    gc = gca;
    hold on
    scatter(valueFix, gc.YLim(1), 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', [1,0,0])
    hold off
    xlabel(p, 'Interpreter', 'none')
    ylabel('cost')
end
i = i + 2;
sp = subplot(nrows, ncols, i); plot(1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan, 1, nan); set(sp, 'Visible', 'off');

legLabel = costFields;
legend(sp, legLabel,'Location', 'west', 'Interpreter', 'none');


% And now displayng mean cost ascribable to nutrient data or size data
plt3 = figure;
plt3.Units = 'inches';
plt3.Position = [0 0 16 8];

nrows = floor(sqrt(npars));
ncols = ceil(npars / nrows);

scalarInd = contains(costFields, Data.scalar.obsInCostFunction);
sizeInd = contains(costFields, Data.size.obsInCostFunction);

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
    cost = squeeze(costComponents(:,ind,:));
    
    costScalar = mean(cost(:,scalarInd), 2);
    costSize = mean(cost(:,sizeInd), 2);
   
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




% Display marginal sensitivities -- plot needs improvement...
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


close all

%% Ensemble outputs

% Generate model output for an ensemble of parameter sets, varying only
% a subest of parameters.

% Within any ensemble, the number of parameters that may be varied should
% probably not exceed 5 or 6, and using only 3 or 4 would be better. This
% is due to the 'curse of dimensionality'... varying lots of parameters
% over a grid and generating model output for each combination, very
% quickly becomes unfeasible as the number of parameters increases...

loadEnsemble = true; % load prior run

switch loadEnsemble
    
    case true
        
        % Load stored results
        fileName = 'paramEnsemble';
        fileName = [fileName, fileSuffix];
        m2 = matfile(fullfile(Directories.resultsDir, fileName), 'Writable', true);
        ParamEnsemble = m2.results;
        vPars = ParamEnsemble.parNames;
        npars = length(vPars);
        parRes = size(ParamEnsemble.value, 1);
        nVals = parRes ^ npars;

    case false

        vPars = {'wPOM1', 'rDON', 'rPON'}; % Parameters to vary
        parNames = FixedParams.tunePars; % Any tuning parameters may be selected, and this should also work for any other model parameters...
        parRes = 15; % Number of values per parameter -- the resolution
        npars = length(vPars);
        nVals = parRes ^ npars; % number of unique combinations -- how many times the model needs evaluated
        parVals = nan(parRes, npars); % all parameter values
        
        for i = 1:npars
            j = strcmp(parNames, vPars{i});
            parVals(:,i) = linspace(boundsLower(j), boundsUpper(j), parRes);
        end
        
        % Store all outputs in ParamEnsemble struct
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
        close(progress)

end


saveOutput = true;

switch saveOutput, case true
    fileName = 'paramEnsemble';
    fileName = [fileName, fileSuffix];
    m2 = matfile(fullfile(Directories.resultsDir, fileName), 'Writable', true);
    m2.results = ParamEnsemble;
end







% Use ensemble model output to assess nutrient component of cost -- using
% various different functions.

% selectFunction = 'meanCDFdist_Hellinger';
% selectFunction = 'smoothCDFdist_Hellinger';
% selectFunction = 'LeastAbsErr_Hellinger';

selectFunction = 'RMS_Hellinger2';


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
close(progress)



% Make surface or contour plots to display cost function values relative to
% two parameters -- slices out of the overall cost surface.

costFields = ParamEnsemble.costFields;

scalarInd = contains(costFields, Data.scalar.obsInCostFunction);
sizeInd = contains(costFields, Data.size.obsInCostFunction);

cost = ParamEnsemble.cost;
costComponents = ParamEnsemble.costComponents;
costTot = cost;
costNut = mean(costComponents(:,scalarInd), 2, 'omitnan'); % careful here... the cost function may average differently...
costSize = mean(costComponents(:,sizeInd), 2, 'omitnan');

newDim = repmat(parRes, [1, npars]); % reshape dimension
allValues = ParamEnsemble.allValues;
allValues = reshape(allValues, [newDim, npars]);
costTot = reshape(costTot, newDim);
costNut = reshape(costNut, newDim);
costSize = reshape(costSize, newDim);



% This plotting code below obviously needs work... but it's useful for a
% detailed look at cost function sensitivities to parameter values, and 
% parameter correlations...

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



