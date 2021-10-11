function [cost, costComponents, costFunctionChoices, costLabel_dataType] = costFunction(varargin)

% Labels for cost function options -- code & comments for each below.
% To use a different cost function simply write its label here then include
% another case below -- it must be fully implementable from Data and
% modData. Note below the requirement of cost function label when fitting
% to full (unbinned) size spectra.
costFunctionChoices = { ...
    'LSS'; ...                                                              % requires integrated/binned size data
    'RMS'; ...                                                              % requires integrated/binned size data
    'Hellinger_groupWaterOrigin'; ...                                       % requires integrated/binned size data
    'IQD_Hellinger_groupWaterOrigin'; ...                                   % requires integrated/binned size data
    'RMS_Hellinger'; ...                                                    % requires integrated/binned size data
    'RMS_Hellinger_ZPratio'; ...                                            % requires integrated/binned size data
    'RMSsmooth_Hellinger'; ...                                              % requires integrated/binned size data
    'RMS_HellingerFullSpectrum'; ...                                        % requires full size spectra
    'RMS_HellingerFullSpectrum_averagedEventsDepths'                        % requires full size spectra
    };

% Specify type of size data used for each cost function choice: either
% 'binnedSizeData' or 'fullSizeSpectra'.
% To do this automatically we require the cost function label to contain
% the string 'FullSpectrum' iff it's designed to fit to the unbinned size
% spectra data.
dataTypes = {'binnedSizeData', 'fullSizeSpectra'};
di = 1 + cellfun(@(z) contains(z, 'FullSpectrum'), costFunctionChoices);
costLabel_dataType = cell(length(costFunctionChoices), 2);
costLabel_dataType(:,1) = costFunctionChoices;
costLabel_dataType(:,2) = arrayfun(@(z) dataTypes(z), di);
costLabel_dataType = cell2table(costLabel_dataType, ... 
    'VariableNames', {'label','dataType'});

% Evaluating cost requires passing name-value pairs for 'label', 'Data' and
% 'modData' as optional varargin values. The 'label' specifies which cost
% function to use; 'Data' are measurements; 'modData' are the modelled
% equivalents.
extractVarargin(varargin)

calculateCost = exist('Data', 'var') && exist('modData', 'var') && exist('label', 'var');
if ~calculateCost
    cost = nan;
    costComponents = nan;
    return
end

label = eval('label');
if ~ismember(label, costFunctionChoices)
    error('Cost function name must match one of the available options in costFunction.m.')
end

% Data types used to fit the model
Vars = Data.scalar.obsInCostFunction;
VarsSize = Data.size.obsInCostFunction;

%% Evaluate cost function with name given by 'label'
switch label
    
    case 'LSS'
        % Least sum of squares method on all data types
        
        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            squaredError = (yobs - ymod) .^ 2;
            costComponents.(varLabel) = sum(squaredError(:)) / numel(squaredError);
        end
        
        % Size spectra data
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            % autotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            yobs = Data.size.dataBinned.scaled_Value(ind);
            if isfield(modData.size, 'scaled_Value')
                ymod = modData.size.scaled_Value(ind,:);
            else
                ymod = modData.size.Value(ind,:);
                ymod = (ymod - mean(ymod)) ./ std(ymod);
            end
            squaredError = (yobs - ymod) .^ 2;
            costComponents.([varLabel '_autotroph']) = sum(squaredError(:)) / numel(squaredError);
            % heterotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
            yobs = Data.size.dataBinned.scaled_Value(ind);
            if isfield(modData.size, 'scaled_Value')
                ymod = modData.size.scaled_Value(ind,:);
            else
                ymod = modData.size.Value(ind,:);
                ymod = (ymod - mean(ymod)) ./ std(ymod);
            end
            squaredError = (yobs - ymod) .^ 2;
            costComponents.([varLabel '_heterotroph']) = sum(squaredError(:)) / numel(squaredError);
        end
        fields = fieldnames(costComponents);
        x = struct2array(costComponents);
        scalarCosts = x(contains(fields, Vars));
        sizeCosts = x(contains(fields, VarsSize));
        cost = mean(scalarCosts) + mean(sizeCosts);
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'RMS'
        % Root mean square for all data types
        
        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            absError = abs(yobs - ymod);
            costComponents.(varLabel) = sum(absError(:)) / numel(absError);
        end
        % Size spectra data
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            % autotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            yobs = Data.size.dataBinned.scaled_Value(ind);
            if isfield(modData.size, 'scaled_Value')
                ymod = modData.size.scaled_Value(ind,:);
            else
                ymod = modData.size.Value(ind,:);
                ymod = (ymod - mean(ymod)) ./ std(ymod);
            end
            absError = abs(yobs - ymod);
            costComponents.([varLabel '_autotroph']) = sum(absError(:)) / numel(absError);
            % heterotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
            yobs = Data.size.dataBinned.scaled_Value(ind);
            if isfield(modData.size, 'scaled_Value')
                ymod = modData.size.scaled_Value(ind,:);
            else
                ymod = modData.size.Value(ind,:);
                ymod = (ymod - mean(ymod)) ./ std(ymod);
            end
            absError = abs(yobs - ymod);
            costComponents.([varLabel '_heterotroph']) = sum(absError(:)) / numel(absError);
%             costComponents.(varLabel) = costComponents.([varLabel '_autotroph']) + costComponents.([varLabel '_heterotroph']);
        end
        fields = fieldnames(costComponents);
        x = struct2array(costComponents);
        scalarCosts = x(contains(fields, Vars));
        sizeCosts = x(contains(fields, VarsSize));
        cost = mean(scalarCosts) + mean(sizeCosts);
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'Hellinger_groupWaterOrigin'
        % To improve the compatability of cost components between different
        % data types use Hellinger distances for scalar (nutrient) and vector
        % (size) data. The scalar data is already transformed to standard
        % normal distribution so using a distribution-based cost metric
        % ought to be workable...
        % In this function a different, but comparable, metric is used for
        % the total abundance info contained in the size data. This is
        % necessary because Hellinger distance is not applicable to these
        % single data points. The metric also returns values in the
        % interval [0,1], with 0 represetning a perfect fit and 1 an very
        % poor fit.        
        % THERE ARE ISSUES WITH THIS FUNCTION: PDFs FOR THE SCALAR DATA ARE
        % CALCULATED NUMERICALLY, WHICH REQUIRES CHOSING THE RESOLUTION OF
        % THE GRID OVER WHICH PDFs ARE GENERATED. THE GRID RESOLUTION 
        % (ARBITRARY CHOICE) AFFECTS THE COST... NOT ACCEPTABLE!
        % A POSSIBLE SOLUTION, KEEPING THE 'SPIRIT' OF THIS METHOD, WOULD
        % BE TO USE THE CDFs OF THE SCALAR DATA RATHER THAN PDFs.
        % INTERPOLATION WOULD BE REQUIRED TO ALIGN THE MODELLED CDFs WITH
        % OBSERVED CDFs, AND THE METRIC WOULD BE SOME AVERAGE BETWEEN-CDF
        % DISTANCE, RATHER THAN HELLINGER DISTANCE, BUT IT WOULD BE A
        % CLEANER (LESS BIASED BY ARBITRARY CHOICE) METHOD...
        
        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);
            m = size(ymod, 2);
            % As the standardised data are approximately normally distributed,
            % we could assume that identically transformed model outputs
            % are also normally distributed, then find their means and
            % variances to calculate Hellinger distances analytically.
            % However, assuming normality of transformed model output is
            % inappropriate because it could be any shape...
            % Thus, derive cdfs and pdfs directly without assuming that
            % transformed model outputs follow normal distributions.
            % To compare observed and modelled distributions, the pdfs will
            % need to be evaluated across identical domains.
            cdf = (1:n)' ./ n;
            yobsc = sort(yobs); % sort observations
            ymodc = sort(ymod);
            % remove any duplicate values (more likely in model output than
            % in data), retaining the largest probability values
            keep = ~[diff(yobsc) == 0; false];
            cdfobs = cdf(keep);
            yobsc = yobsc(keep);
            % store cdfmod as cell-array because vector lengths may differ
            % after removing duplicates
            keep = ~[diff(ymodc) == 0; false(1, m)];
            cdfmod = cell(1, m);
            ymodc_ = cell(1, m);
            for ij = 1:m
                cdfmod{:,ij} = cdf(keep(:,ij));
                ymodc_{:,ij} = ymodc(keep(:,ij),ij);
            end
            % define regular grid across measurement space -- define number
            % of grid nodes, ng, relative to number of observations, n.
            vrange = [min([yobsc(:); ymodc(:)]), max([yobsc(:); ymodc(:)])];
            % choice of grid resolution influences cost values... is there a better way to do this, or a principled way to choose grid resolution?
            nr = 1; % grid resolution (0 < nr <= 1)
            ng = ceil(nr * n);
            grid = nan(ng, 1);
            grid(2:end) = linspace(vrange(1), vrange(2), ng-1)';
            sp = diff(grid(2:3));
            grid(1) = grid(2) - sp;
            % interpolate cdfs and derive pdfs
            cdfobsi = interp1(yobsc, cdfobs, grid, 'linear', 'extrap');
            cdfobsi(cdfobsi < min(cdfobs)) = 0;
            cdfobsi(cdfobsi > 1) = 1;
            pdfobsi = diff(cdfobsi);
            cdfmodi = nan(ng, m);
            pdfmodi = nan(ng-1, m);
            for ij = 1:m
                cdfmodi(:,ij) = interp1(ymodc_{ij}, cdfmod{ij}, grid, 'linear', 'extrap');
                cdfmodi(cdfmodi(:,ij) < min(cdfmod{ij}), ij) = 0;
                cdfmodi(cdfmodi(:,ij) > 1, ij) = 1;
                pdfmodi(:,ij) = diff(cdfmodi(:,ij));
            end
            hellingerDistance = (1 - sum((pdfobsi .* pdfmodi) .^ 0.5)) .^ 0.5;
            costComponents.(varLabel) = mean(hellingerDistance); % average over trajectory selections
        end
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Average the cost across data-types (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'IQD_Hellinger_groupWaterOrigin'
        
        % This function addresses the grid resolution problem of case
        % Hellinger_groupWaterOrigin by using a CDF-based metric for the
        % scalar data -- the integrated quadratic distance between observed
        % and modelled CDFs.
        % Using differences between CDFs it is difficult to define a metric
        % bounded in [0,1] that is sufficiently robust to be a useful
        % measure of cost. When the model values are outside the range of
        % the data, differences between CDFs may be insensitive to changes
        % in parameter values => the cost may be insensitive. This is an
        % issue only for the scalar (nutrient) data because the domain of 
        % the observed and modelled distributions are not guarenteed to
        % overlap -- unlike the vector (size) data which has consistent
        % domain (size classes) regardless of numeric values.

        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);
            m = size(ymod, 2);
            yobsc = sort(yobs); % sort observations
            ymodc = sort(ymod);
            cdf = (1:n)' ./ n;
            % remove any duplicate values (more likely in model output than
            % in data), retaining the largest probability values
            keep = ~[diff(yobsc) == 0; false];
            cdfobs = cdf(keep);
            yobsc = yobsc(keep);
            % store cdfmod as cell-array because vector lengths may differ
            % after removing duplicates
            keep = ~[diff(ymodc) == 0; false(1, m)];
            cdfmod = cell(1, m);
            ymodc_ = cell(1, m);
            for ij = 1:m
                cdfmod{:,ij} = cdf(keep(:,ij));
                ymodc_{:,ij} = ymodc(keep(:,ij),ij);
            end
            
            % interpolate modelled CDFs over query points defined by
            % the data
            cdfmodi = cell(1,m);
            for ij =1:m
                cdfmodi{ij} = interp1(ymodc_{ij}, cdfmod{ij}, yobsc, 'linear', 'extrap');
                cdfmodi{ij}(cdfmodi{ij} < 0) = 0;
                cdfmodi{ij}(cdfmodi{ij} > 1) = 1;
            end
            
            % distances between each data CDF value and equivalent
            % modelled CDF value
            cdfDist = cell(1,m);
            avDist = nan(1,m);
            for ij = 1:m
                cdfDist{ij} = abs(cdfobs - cdfmodi{ij});
                avDist(ij) = trapz(cdfDist{ij}) / length(cdfDist{ij});
%                 avDist(ij) = mean(cdfDist{ij});
            end
            costComponents.(varLabel) = mean(avDist); % average over trajectory selections
        end
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end

        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Average the cost across data-types (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'RMS_Hellinger'
        
        % Root-mean-squared (absolute) errors for scalar/nutrient
        % data.  Cost is independent of arbitrary assumptions and it
        % remains sensitive to parameters for any model outputs (ie, 
        % addresses issues with Hellinger_groupWaterOrigin and 
        % IQD_Hellinger_groupWaterOrigin). However, unlike the CDF-based
        % functions above, outputs are not bounded in [0,1], but if the
        % data range is used as a scaling term then a good model fit to
        % data produces cost values of the same order as the Hellinger
        % distances ([0,1]) used for the relative abundance-at-size data.
        
        % Scalar data
        plotFit = false;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);            
            % Sort observations
            [yobs_sort, o] = sort(yobs);
            CDFobs = [yobs_sort, ...
                (1:n)' ./ n]; % Empirical CDF of standardised data
            % Reorder modelled values to match sorted data
            ymod_sort = ymod(o,:);
            CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
            switch plotFit, case true
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                hold on
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
            end
            %~
            % This has been modified from a root-mean-square (RMS) cost.
            % True RMS function code is commented.
            %~
%             err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
            err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
%             avFun = @mean;
            avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
%             RMS = (avFun(err2)) .^ 0.5; % root-mean-square
            RMS = avFun(err); % average absolute error
            RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
            
            costComponents.(varLabel) = mean(RMS); % average over trajectory selections
            
        end
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
            end
        end
        
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Treat PON & POC data as a single type, POM, thereby
        % downweighting their combined cost contribution
        POMi = ismember(Vars, {'PON','POC'});
        costPOM = mean(costNutrient(POMi));
        costNutrient(POMi) = [];
        costNutrient = [costNutrient, costPOM];
        
        % Average the cost across data-types (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components
        %             cost = cost(2); % try fitting only to the size data
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    case 'RMS_Hellinger_ZPratio'
        
        % As above, but for the total abundance in size data fit the
        % ratio that's zooplankton. This means that the chlorophyll
        % data should inform the phytoplankton abundance while the
        % zooplankton abundance is constrained relative to
        % phytoplankton.
        
        % Scalar data
        plotFit = false;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);
            % Sort observations
            [yobs_sort, o] = sort(yobs);
            CDFobs = [yobs_sort, ...
                (1:n)' ./ n]; % Empirical CDF of standardised data            
            % Reorder modelled values to match sorted data
            ymod_sort = ymod(o,:);
            CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
            switch plotFit, case true
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                hold on
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
            end
            %~
            % This has been modified from a root-mean-square (RMS) cost.
            % True RMS function code is commented.
            %~
%             err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
            err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
%             avFun = @mean;
            avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
%             RMS = (avFun(err2)) .^ 0.5; % root-mean-square
            RMS = avFun(err); % average absolute error
            RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
            
            costComponents.(varLabel) = avFun(RMS); % average over trajectory selections
            
        end
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTotP = sum(yobs);
                ymodTotP = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTotP;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTotP;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTotZ = sum(yobs);
                ymodTotZ = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTotZ;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTotZ;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                
                % Cost value for heterotroph total biovolume
                yobsTot = yobsTotP + yobsTotZ;
                ymodTot = ymodTotP + ymodTotZ;
                cost_totZ = abs((yobsTotZ ./ yobsTot) - (ymodTotZ ./ ymodTot)); % absolute difference of probabilities (bounded in (0,1))
                
                costComponents.(label_autoTot) = nan;
                costComponents.(label_heteroTot) = mean(cost_totZ);
                
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 1;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
%                 cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
%                 costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,1,i) = cra;
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Treat PON & POC data as a single type, POM, thereby
        % downweighting their combined cost contribution
        POMi = ismember(Vars, {'PON','POC'});
        costPOM = mean(costNutrient(POMi));
        costNutrient(POMi) = [];
        costNutrient = [costNutrient, costPOM];
        
        % Average the cost for each data type (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components
        
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'RMSsmooth_Hellinger'
        % Same as RMS_Hellinger but modelled output are smoothed to reduce
        % overfitting outliers.

        
        % Scalar data
        smoothFactor = 0.35; % smoothFactor multiplies number of data points to give loess smoothing parameter
        plotFit = false;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);            
            % Sort observations
            [yobs_sort, o] = sort(yobs);
            CDFobs = [yobs_sort, ...
                (1:n)' ./ n]; % Empirical CDF of standardised data
            % Reorder modelled values to match sorted data
            ymod_sort = ymod(o,:);
            CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
            switch plotFit, case true
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                hold on
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
            end
            % Smooth model output
            smoothF = smoothFactor * n;
            xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothF, 'loess');
            %~
            % This has been modified from a root-mean-square (RMS) cost.
            % True RMS function code is commented.
            %~
%             err2 = (CDFobs(:,1) - xx) .^ 2; % squared error
            err = abs(CDFobs(:,1) - xx); % absolute error
%             avFun = @mean;
            avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
%             RMS = (avFun(err2)) .^ 0.5; % root-mean-square
            RMS = avFun(err); % average absolute error
            RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
            
            costComponents.(varLabel) = mean(RMS); % average over trajectory selections
            
        end
        
        % Vector (size) data
        groupedByWaterOrigin = isfield(Data.size.dataBinned, 'waterMass');
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.dataBinned.waterMass);
            else
                waterMasses = {[]};
            end
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0;
                if groupedByWaterOrigin
                    ind1 = ind0 & strcmp(Data.size.dataBinned.waterMass, wm);
                end
                if groupedByWaterOrigin
                    label_autoRel = [varLabel '_' wm '_autotroph_Rel'];
                    label_autoTot = [varLabel '_' wm '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_' wm '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_' wm '_heterotroph_Tot'];
                else
                    label_autoRel = [varLabel '_autotroph_Rel'];
                    label_autoTot = [varLabel '_autotroph_Tot'];
                    label_heteroRel = [varLabel '_heterotroph_Rel'];
                    label_heteroTot = [varLabel '_heterotroph_Tot'];
                end
                
                % autotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_autoRel) = mean(hellingerDistance); % average over trajectory selections
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                
                % heterotrophs
                ind = ind1 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
                yobs = Data.size.dataBinned.Value(ind);
                ymod = modData.size.Value(ind,:);
                nsize = size(ymod, 1);
                ymod = reshape(ymod(~isnan(ymod)), nsize, []); % remove NaNs ymod contains when fitting to Arctic AND Atlantic data
                nsample = size(ymod, 2);
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                % Cost metric for relative abundance (vector)
                cdfobs = cumsum(yobs) ./ yobsTot;
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ ymodTot;
                pdfmod = diff([zeros(1, nsample); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.(label_heteroRel) = mean(hellingerDistance);
                % Cost metric for total abundance (scalar)
                u = abs(log(ymodTot / yobsTot));
                u = exp(-a .* u);
                z = (1 - u) ./ (1 + u);
                costComponents.(label_heteroTot) = mean(z); % average over trajectory selections
            end
        end
        
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Treat PON & POC data as a single type, POM, thereby
        % downweighting their combined cost contribution
        POMi = ismember(Vars, {'PON','POC'});
        costPOM = mean(costNutrient(POMi));
        costNutrient(POMi) = [];
        costNutrient = [costNutrient, costPOM];
        
        % Average the cost across data-types (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components
        %             cost = cost(2); % try fitting only to the size data
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    case 'RMS_HellingerFullSpectrum'
        
        % Same as case RMS_Hellinger but use model output representing the
        % full size spectra rather than binned (integrated) data.
        % Build the size data cost section to separately fit all size
        % spectra -- keep the returned cost components as autotrophs &
        % heterotrophs in Arctic and/or Atlantic waters => average the
        % costs over sampling events and depths within this function.
        
        % Scalar data
        plotFit = false;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);
            % Sort observations
            [yobs_sort, o] = sort(yobs);
            CDFobs = [yobs_sort, ...
                (1:n)' ./ n]; % Empirical CDF of standardised data            
            % Reorder modelled values to match sorted data
            ymod_sort = ymod(o,:);
            CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
            switch plotFit, case true
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                hold on
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
            end
            %~
            % This has been modified from a root-mean-square (RMS) cost.
            % True RMS function code is commented.
            %~
%             err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
            err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
%             avFun = @mean;
            avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
%             RMS = (avFun(err2)) .^ 0.5; % root-mean-square
            RMS = avFun(err); % average absolute error
            RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
            
            costComponents.(varLabel) = avFun(RMS); % average over trajectory selections
            
        end
        
        % Vector (size) data
        trophicLevels = unique(Data.size.trophicLevel);
        groupedByWaterOrigin = isfield(Data.size, 'waterMass');
        if groupedByWaterOrigin
            waterMasses = unique(Data.size.waterMass);
        else
            waterMasses = {[]};
        end
        % VarsSize labels match binned data => redefine for full size spectra
        switch Data.size.obsInCostFunction{1}
            case 'BioVol'; VarsSize = {'BioVolDensity'};
            case 'CellConc'; VarsSize = {'cellDensity'};
        end
        allEvents = unique(Data.size.Event, 'stable');
        nVars = length(VarsSize);
        nWaterMasses = length(waterMasses);
        nTrophicLevels = length(trophicLevels);
        HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
        totAbnMetric = nan(nVars, nWaterMasses, nTrophicLevels);
        a = log(3) / log(2); % steepness of cost metric for total abundance
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind0 = ismember(Data.size.Event, allEvents(strcmp(Data.size.waterMass, wm)));
                Events = unique(Data.size.Event(ind0));
                counter = 0;
                for ie = 1:length(Events)
                    event = Events(ie);
                    ind1 = ind0 & Data.size.Event == event;
                    depths = unique(Data.size.Depth(ind1));
                    for id = 1:length(depths)
                        counter = counter + 1;
                        depth = depths(id);
                        ind2 = ind1 & Data.size.Depth == depth;
                        for it = 1:length(trophicLevels)
                            trophicLevel = trophicLevels{it};
                            ind3 = ind2 & strcmp(Data.size.trophicLevel, trophicLevel);
                            yobs = Data.size.(varLabel)(ind3);
                            ymod = modData.size.(varLabel)(ind3);
                            % Hellinger distance for relative abundance-at-size                            yobsRel = yobs ./ sum(yobs);
                            yobsRel = yobs ./ sum(yobs);
                            ymodRel = ymod ./ sum(ymod);
                            HellingerDistance = (1 - sum((yobsRel .* ymodRel) .^ 0.5)) .^ 0.5;
                            HellingerDistances(i,w,it,counter) = HellingerDistance;
                            % Metric for total abundance
                            yobsTot = trapz(yobs);
                            ymodTot = trapz(ymod);
                            u = abs(log(ymodTot / yobsTot));
                            u = exp(-a .* u);
                            z = (1 - u) ./ (1 + u);
                            totAbnMetric(i,w,it,counter) = z;
                        end
                    end
                end
            end
        end
        
        % Average Hellinger distances over sample events & depths
        HellingerDistances_ = mean(HellingerDistances, ndims(HellingerDistances));
        totAbnMetric_ = mean(totAbnMetric, ndims(totAbnMetric));
        % Store in costComponents
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for it = 1:length(trophicLevels)
                trophicLevel = trophicLevels{it};
                if groupedByWaterOrigin
                    for w = 1:length(waterMasses)
                        waterMass = waterMasses{w};
                        label = [varLabel '_' waterMass '_' trophicLevel];
                        costComponents.([label '_Rel']) = HellingerDistances_(i,w,it);
                        costComponents.([label '_Tot']) = totAbnMetric_(i,w,it);
                    end
                else
                    label = [varLabel '_' trophicLevel];
                    costComponents.([label '_Rel']) = HellingerDistances_(i,1,it);
                    costComponents.([label '_Tot']) = totAbnMetric_(i,1,it);
                end
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Treat PON & POC data as a single type, POM, thereby
        % downweighting their combined cost contribution
        POMi = ismember(Vars, {'PON','POC'});
        costPOM = mean(costNutrient(POMi));
        costNutrient(POMi) = [];
        costNutrient = [costNutrient, costPOM];
        
        % Average the cost across data-types (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    case 'RMS_HellingerFullSpectrum_averagedEventsDepths'
        
        % Same as RMS_HellingerFullSpectrum_averagedEventsDepths but
        % average over all events and depths before comparing model to data
        
        % Scalar data
        plotFit = false;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);
            % Sort observations
            [yobs_sort, o] = sort(yobs);
            CDFobs = [yobs_sort, ...
                (1:n)' ./ n]; % Empirical CDF of standardised data            
            % Reorder modelled values to match sorted data
            ymod_sort = ymod(o,:);
            CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
            switch plotFit, case true
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                hold on
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
            end
            %~
            % This has been modified from a root-mean-square (RMS) cost.
            % True RMS function code is commented.
            %~
%             err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
            err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
%             avFun = @mean;
            avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
%             RMS = (avFun(err2)) .^ 0.5; % root-mean-square
            RMS = avFun(err); % average absolute error
            RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
            
            costComponents.(varLabel) = avFun(RMS); % average over trajectory selections
            
        end
        
        % Vector (size) data
        trophicLevels = unique(Data.size.trophicLevel);
        groupedByWaterOrigin = isfield(Data.size, 'waterMass');
        if groupedByWaterOrigin
            waterMasses = unique(Data.size.waterMass);
        else
            waterMasses = {[]};
        end
        % VarsSize labels match binned data => redefine for full size spectra
        switch Data.size.obsInCostFunction{1}
            case 'BioVol'; VarsSize = {'BioVolDensity'};
            case 'CellConc'; VarsSize = {'cellDensity'};
        end
        allEvents = unique(Data.size.Event, 'stable');
        nVars = length(VarsSize);
        nWaterMasses = length(waterMasses);
        nTrophicLevels = length(trophicLevels);
        nESD = length(unique(Data.size.ESD));
        a = log(3) / log(2); % steepness of cost metric for total abundance
        % Extract modelled output and data into arrays with sample
        % event and depth in the trailing dimension
        yobsRel = nan(nVars, nWaterMasses, nTrophicLevels, nESD);
        ymodRel = nan(nVars, nWaterMasses, nTrophicLevels, nESD);
        yobsTot = nan(nVars, nWaterMasses, nTrophicLevels);
        ymodTot = nan(nVars, nWaterMasses, nTrophicLevels);
        for i = 1:nVars
            varLabel = VarsSize{i};
            for w = 1:nWaterMasses
                wm = waterMasses{w};
                ind0 = ismember(Data.size.Event, allEvents(strcmp(Data.size.waterMass, wm)));
                Events = unique(Data.size.Event(ind0));
                nEvents = length(Events);
                for it = 1:nTrophicLevels
                    trophicLevel = trophicLevels{it};
                    ind1 = ind0 & strcmp(Data.size.trophicLevel, trophicLevel);
                    counter = 0;
                    for ie = 1:nEvents
                        event = Events(ie);
                        ind2 = ind1 & Data.size.Event == event;
                        Depths = unique(Data.size.Depth(ind2));
                        nDepths = length(Depths);
                        for id = 1:nDepths
                            counter = counter + 1;
                            Depth = Depths(id);
                            ind3 = ind2 & Data.size.Depth == Depth;
                            yobs = Data.size.(varLabel)(ind3);
                            ymod = modData.size.(varLabel)(ind3);
                            yobsRel(i,w,it,:,counter) = yobs ./ sum(yobs);
                            ymodRel(i,w,it,:,counter) = ymod ./ sum(ymod);
                            yobsTot(i,w,it,counter) = trapz(yobs);
                            ymodTot(i,w,it,counter) = trapz(ymod);
                        end
                    end
                end
            end
        end
        % Average over sample events and depths
        if size(yobsRel, ndims(yobsRel)) ~= nESD
            yobsRel = mean(yobsRel, ndims(yobsRel));
            ymodRel = mean(ymodRel, ndims(ymodRel));
            yobsTot = mean(yobsTot, ndims(yobsTot));
            ymodTot = mean(ymodTot, ndims(ymodTot));
        end
        
        % Find Hellinger distances betwen model output and data for
        % these averaged spectra, for all variables, water masses and
        % trophic levels.
        HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
        totAbnMetric = nan(nVars, nWaterMasses, nTrophicLevels);
        for i = 1:nVars
            for w = 1:nWaterMasses
                for it = 1:nTrophicLevels
                    yobs_ = yobsRel(i,w,it,:);
                    ymod_ = ymodRel(i,w,it,:);
                    HellingerDistance = (1 - sum((yobs_ .* ymod_) .^ 0.5)) .^ 0.5;
                    HellingerDistances(i,w,it) = HellingerDistance;
                    yobsTot_ = yobsTot(i,w,it);
                    ymodTot_ = ymodTot(i,w,it);
                    u = abs(log(ymodTot_ / yobsTot_));
                    u = exp(-a .* u);
                    z = (1 - u) ./ (1 + u);
                    totAbnMetric(i,w,it) = z;
                end
            end
        end
                
        % Store in costComponents
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            for it = 1:length(trophicLevels)
                trophicLevel = trophicLevels{it};
                if groupedByWaterOrigin
                    for w = 1:length(waterMasses)
                        waterMass = waterMasses{w};
                        label = [varLabel '_' waterMass '_' trophicLevel];
                        costComponents.([label '_Rel']) = HellingerDistances(i,w,it);
                        costComponents.([label '_Tot']) = totAbnMetric(i,w,it);
                    end
                else
                    label = [varLabel '_' trophicLevel];
                    costComponents.([label '_Rel']) = HellingerDistances(i,1,it);
                    costComponents.([label '_Tot']) = totAbnMetric(i,1,it);
                end
            end
        end
        
        % Within each size data group, weight relative abundance-at-size
        % relative to total abundance.
        weight_relVsTot = 3;
        nDataTypes = length(VarsSize);
        costSize = zeros(nDataTypes, 2, length(waterMasses)); % store weighted costs for all size data groups (dimension = data type, trophic level, water mass)
        for j = 1:nDataTypes
            varLabel = VarsSize{j};
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cLabel = [varLabel '_' wm '_'];
                else
                    cLabel = [varLabel '_'];
                end
                cta = costComponents.([cLabel 'autotroph_Tot']);
                cra = costComponents.([cLabel 'autotroph_Rel']);
                cth = costComponents.([cLabel 'heterotroph_Tot']);
                crh = costComponents.([cLabel 'heterotroph_Rel']);
                costSize(j,1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(j,2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            end
        end
        
        costNutrient = zeros(1,length(Vars));
        for i = 1:length(Vars)
            costNutrient(i) = costComponents.(Vars{i});
        end
        
        % Treat PON & POC data as a single type, POM, thereby
        % downweighting their combined cost contribution
        POMi = ismember(Vars, {'PON','POC'});
        costPOM = mean(costNutrient(POMi));
        costNutrient(POMi) = [];
        costNutrient = [costNutrient, costPOM];
        
        % Average the cost across data-types (nutrient & size)
        cost = [mean(costNutrient), mean(costSize(:))];
        % Assign size vs nutrient weighting
        weight_sizeVsNutrient = 1; % weighting factor of size vs nutrient data
        weights = 2 .* [1, weight_sizeVsNutrient] ./ (weight_sizeVsNutrient+1);
        cost = weights .* cost;
        
        cost = mean(cost); % finally, average cost over nutrient and size data components

end


%%

if ~exist('cost', 'var')
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if ~exist('costComponents', 'var')
    costComponents = nan;
end
