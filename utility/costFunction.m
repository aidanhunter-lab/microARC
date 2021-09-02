function [cost, costComponents, costFunctionChoices] = costFunction(varargin)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Labels for cost function options -- code for each is below.
% To use a different cost function simply write its label here then include
% another case below -- it must be fully code-able from Data and modData
costFunctionChoices = { ...
    'LSS', ...
    'RMS', ...
    'Hellinger_groupWaterOrigin', ...
    'Hellinger2_groupWaterOrigin', ...
    'Hellinger3_groupWaterOrigin', ...
    'IQD_Hellinger_groupWaterOrigin', ...
    'meanCDFdist_Hellinger', ...
    'smoothCDFdist_Hellinger', ...
    'LeastAbsErr_Hellinger', ...
    'RMSsmooth_Hellinger', ...
    'RMS_Hellinger', ...
    'RMS_Hellinger2', ...
    'meanCDFdist_HellingerFullSpectrum', ...
    'meanCDFdist_HellingerFullSpectrum_averagedEventsDepths'
    };

% Evaluating cost requires passing name-value pairs for 'label', 'Data' and
% 'modData' as optional varargin values. The 'label' specifies which cost
% function to use; 'Data' are measurements; 'modData' are the modelled
% equivalents.
extractVarargin(varargin)
calculateCost = exist('Data', 'var') & exist('modData', 'var') & exist('label', 'var');

if ~calculateCost
    
    cost = nan;
    costComponents = nan;
    
else
    
    if ~ismember(label, costFunctionChoices)
        error('Cost function name must match one of the available options in costFunction.m.')
    end
    
    % Data types used to fit the model
    Vars = Data.scalar.obsInCostFunction;
    VarsSize = Data.size.obsInCostFunction;

    switch label
        
        case 'LSS' % Least sum of squares
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
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
            
        case 'RMS' % Root mean square
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                absError = sqrt((yobs - ymod) .^ 2);
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
                absError = sqrt((yobs - ymod) .^ 2);
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
                absError = sqrt((yobs - ymod) .^ 2);
                costComponents.([varLabel '_heterotroph']) = sum(absError(:)) / numel(absError);
                costComponents.(varLabel) = costComponents.([varLabel '_autotroph']) + costComponents.([varLabel '_heterotroph']);
            end
            fields = fieldnames(costComponents);
            x = struct2array(costComponents);
            scalarCosts = x(contains(fields, Vars));
            sizeCosts = x(contains(fields, VarsSize));
            cost = mean(scalarCosts) + mean(sizeCosts);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            

        case 'Hellinger2_groupWaterOrigin'
            % Calculate non-parameteric Hellinger distances for the scalar
            % (nutrient) data and for size spectra vectors of relative
            % abundance.
            % Use a different, but comparable, metric for the total abundance
            % data contained in the size spectra. This is necessary because
            % these data do not form distributions so Hellinger distance is not
            % applicable. Like the Hellinger distance, the metric takes values
            % in the in the interval [0,1], with 0 represetning a perfect fit
            % and 1 being an extremely poor fit.
            % This method makes zero assumptions about variability -- unlike
            % the other methods which use variability across particle
            % trajectories to create distributions. The Hellinger distances are
            % calculated numerically instead of using formulas based on
            % assumptions of Gaussian distributions. The scalar data are
            % continuous and approximately Gaussian whereas the size spectra
            % data are discrete, so Hellinger distance calculations are
            % different for each...
            
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
                nr = 1; % choice of grid resolution influences cost values... is there a better way to do this, or a principled way to choose grid resolution?
%                 ng = 9; 
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
            weight_relVsTot = 3; % weighting factor of relative vs total abundance
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
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
            
        case 'Hellinger3_groupWaterOrigin'
            
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
                
                % interpolate data-CDF over modelled query points
                cdfobsi = cell(1,m);
                for ij = 1:m
                    cdfobsi{ij} = interp1(yobsc, cdfobs, ymodc_{ij}, 'linear', 'extrap');
                    cdfobsi{ij}(cdfobsi{ij} < 0) = 0;
                    cdfobsi{ij}(cdfobsi{ij} > 1) = 1;
                end
                % distances between each modelled CDF value and equivalent
                % data-CDF value, and use average distance as the metric
                cdfDist = cell(1,m);
                avDist = nan(1,m);
                for ij = 1:m
                    cdfDist{ij} = abs(cdfobsi{ij} - cdfmod{ij});
                    avDist(ij) = mean(cdfDist{ij});
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
            weight_relVsTot = 3; % weighting factor of relative vs total abundance
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
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
                % modelled CDF value -- use average distance as the metric
                cdfDist = cell(1,m);
                avDist = nan(1,m);
                for ij = 1:m
                    cdfDist{ij} = (cdfobs - cdfmodi{ij}) .^ 2;
                    avDist(ij) = mean(cdfDist{ij}) .^ 0.5;
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
            weight_relVsTot = 3; % weighting factor of relative vs total abundance
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
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

        case 'meanCDFdist_Hellinger'
            
            % This method uses [0,1] metrics for all data types.
            % The scalar data are fitted by reference to their CDFs in a
            % manner that explicitly accounts for each data point rather
            % than fitting to a distribution.
            % The relative abundance-at-size data are fitted with Hellinger
            % distances.
            % The total abundance from the size data are fitted with my
            % custom metric (relies on a ratio rather like Shannon's entropy)
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
%                 m = size(ymod, 2);
                cdf = (1:n)' ./ n;
                [yobs_sort, o] = sort(yobs); % sort observations
                % put modelled output in same order -- these values will
                % almost certainly NOT be in ascending order...
                ymod_sort = ymod(o,:);
%                 figure(1)
%                 plot(yobs_sort, cdf)
%                 hold on
%                 plot(ymod_sort, cdf)
%                 hold off
                % find distances between each modelled data point and the
                % empirical data CDF
                ymodi = interp1(yobs_sort, cdf, ymod_sort, 'linear', 'extrap');
                ymodi = min(1, max(0, ymodi));
                CDFdist = abs(cdf - ymodi);
                avDist = mean(CDFdist); % average over data points
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
            
%             % Within each size data group, weight relative abundance-at-size
%             % relative to total abundance.
%             weight_relVsTot = 3; % weighting factor of relative vs total abundance
%             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
%             for i = 1:length(waterMasses)
%                 wm = waterMasses{i};
%                 if groupedByWaterOrigin
%                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
%                 else
%                     cta = costComponents.([varLabel '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_heterotroph_Rel']);
%                 end
%                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
%             end
                        
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
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
%             cost = cost(2); % try fitting only to the size data
            

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'smoothCDFdist_Hellinger'
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                
%                 % Smooth the CDF
%                 smoothFactor = 0.15;
%                 CDFobs(:,2) = smooth(CDFobs(:,1), CDFobs(:,2), ... 
%                     smoothFactor * n, 'loess');
%                 % Tidy up boundaries
%                 ind = find(CDFobs(:,2) >= 1, 1);
%                 CDFobs(ind:end,2) = 1;
%                 ind = find(CDFobs(:,2) == min(CDFobs(:,2)));
%                 CDFobs(1:ind,2) = min(CDFobs(:,2));
%                 CDFobs(:,2) = max(0, CDFobs(:,2));
%                 plot(CDFobs(:,1), CDFobs(:,2)) % smoothed CDF of standardised data
                
%                 PDFobs = [CDFobs(2:end,1), diff(CDFobs(:,2))];
%                 plot(PDFobs(:,1), PDFobs(:,2))
                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))                
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                % Variability in model outputs may result in a
                % noisy/non-smooth cost function.
                % Smoothing model outputs before comparing to data may help
                % by improving the signal:noise ratio... fitting to pattern
                % in data rather than each individual data point.
                
                % Try fitting by aligning the CDFs of data and model
                
                CDFmod_ = [sort(CDFmod(:,1)), (1:n)' ./ n];
                ind = diff(CDFmod_(:,1)) == 0;
                CDFmod_r = CDFmod_;
                CDFmod_r(ind,:) = [];
                
%                 plot(CDFmod_r(:,1), CDFmod_r(:,2))

                % Evaluate model CDF value for each data point
                matchCDF = interp1(CDFmod_r(:,1), CDFmod_r(:,2), CDFobs(:,1), ...
                    'linear', 'extrap');
                matchCDF = min(1, max(0, matchCDF));
                % Difference between observed and modelled CDFs
                CDFdist = abs(CDFobs(:,2) - matchCDF);
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2))
%                 hold on
%                 plot(CDFmod_(:,1), CDFmod_(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFobs(ij,1)], ...
%                         [CDFobs(ij,2), matchCDF(ij)], 'Color', [0.85, 0.85, 0.85])
%                 end
%                 
                % Average over data points and trajectory selections
                avFun = @mean;
%                 avFun = @geomean;
                avDist = avFun(CDFdist);  % average over data points
                costComponents.(varLabel) = avFun(avDist); % average over trajectory selections
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
            
%             % Within each size data group, weight relative abundance-at-size
%             % relative to total abundance.
%             weight_relVsTot = 3; % weighting factor of relative vs total abundance
%             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
%             for i = 1:length(waterMasses)
%                 wm = waterMasses{i};
%                 if groupedByWaterOrigin
%                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
%                 else
%                     cta = costComponents.([varLabel '_autotroph_Tot']);
%                     cra = costComponents.([varLabel '_autotroph_Rel']);
%                     cth = costComponents.([varLabel '_heterotroph_Tot']);
%                     crh = costComponents.([varLabel '_heterotroph_Rel']);
%                 end
%                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
%                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
%             end
                        
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
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
%             cost = cost(2); % try fitting only to the size data


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'smooth2CDFdist_Hellinger'
            
            % There are problems emerging from using CDF values to
            % calculate the scalar cost values. When modelled values are
            % either much larger or smaller than the data range, the cost
            % becomes insensitive... Probably better to use differences
            % between data and model, rather than converting to CDF
            % values... even though it might be difficult to scale for
            % compatibility with the size-data Hellinger distances.
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
                figure
                plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
                hold on
                
                scatter(CDFmod(:,1), CDFmod(:,2))
                for ij = 1:n
                    plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
                end
                
                % Variability in model outputs may result in a
                % noisy/non-smooth cost function.
                % Smoothing model outputs before comparing to data may help
                % by improving the signal:noise ratio... fitting to pattern
                % in data rather than each individual data point.
                
                smoothFactor = 0.35 * n;
                xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothFactor, 'loess');                
                
                % For each data point, find the smoothed model equivalent,
                % and evaluate the data CDF at those points.
                matchCDF = interp1(CDFobs(:,1), CDFobs(:,2), xx, 'linear', 'extrap');
                matchCDF = min(1, max(0, matchCDF));
                % Difference between observed and modelled CDFs
                CDFdist = abs(CDFobs(:,2) - matchCDF);
                                
                
                figure
                plot(CDFobs(:,1), CDFobs(:,2))
                hold on
                plot(xx, CDFmod(:,2))
                for ij = 1:n
                    plot([xx(ij), xx(ij)], ...
                        [CDFmod(ij,2), matchCDF(ij)], 'Color', [0.85, 0.85, 0.85])
                end
                
                % Average over data points and trajectory selections
%                 avFun = @mean;
                avFun = @geomean;
                avDist = avFun(CDFdist);  % average over data points
                costComponents.(varLabel) = avFun(avDist); % average over trajectory selections
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
            
            %             % Within each size data group, weight relative abundance-at-size
            %             % relative to total abundance.
            %             weight_relVsTot = 3; % weighting factor of relative vs total abundance
            %             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            %             for i = 1:length(waterMasses)
            %                 wm = waterMasses{i};
            %                 if groupedByWaterOrigin
            %                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
            %                 else
            %                     cta = costComponents.([varLabel '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_heterotroph_Rel']);
            %                 end
            %                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
            %                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            %             end
            
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
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
            %             cost = cost(2); % try fitting only to the size data


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'RMS_Hellinger'
            
            % Use root-mean-squared (absolute) errors for scalar/nutrient
            % data. Fit to smoothed RMS values for robustness against
            % overfitting to data points the model cannot replicate
            % (usually the deepest samples).
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                
%                 plot(xx, CDFmod(:,2))
                
                
%                 err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
                err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
                
%                 avFun = @mean;
                avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
                
%                 RMS = (avFun(err2)) .^ 0.5; % average absolute error
                RMS = avFun(err); % average absolute error
                RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
                
                costComponents.(varLabel) = avFun(RMS); % average over trajectory selections

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
            weight_relVsTot = 3; % weighting factor of relative vs total abundance (relative abundance assumed more reliable)
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
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

        case 'RMS_Hellinger2'
            
            % As above, but for the total abundance in size data fit the
            % ratio that's zooplankton. This means that the chlorophyll
            % data should inform the phytoplankton abundance while the
            % zooplankton abundance is constrained relative to
            % phytoplankton.
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                
%                 plot(xx, CDFmod(:,2))
                
                
%                 err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
                err = abs(CDFobs(:,1) - CDFmod(:,1)); % absolute error
                
%                 avFun = @mean;
                avFun = @geomean; % Errors have skewed distributions => use geometric mean for robustness against overfitting data points that the model cannot reproduce
                
%                 RMS = (avFun(err2)) .^ 0.5; % average absolute error
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
%                     % Cost metric for total abundance (scalar)
%                     u = abs(log(ymodTot / yobsTot));
%                     u = exp(-a .* u);
%                     z = (1 - u) ./ (1 + u);
%                     costComponents.(label_autoTot) = mean(z); % average over trajectory selections
                    
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
%                     % Cost metric for total abundance (scalar)
%                     u = abs(log(ymodTot / yobsTot));
%                     u = exp(-a .* u);
%                     z = (1 - u) ./ (1 + u);
%                     costComponents.(label_heteroTot) = mean(z); % average over trajectory selections

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
%             weight_relVsTot = 3; % weighting factor of relative vs total abundance (relative abundance assumed more reliable)
            weight_relVsTot = 1;
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
%                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
%                     cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
%                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(1,i) = cra;
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
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

        case 'RMSsmooth_Hellinger'
            
            % Use root-mean-squared (absolute) errors for scalar/nutrient
            % data. Fit to smoothed RMS values for robustness against
            % overfitting to data points the model cannot replicate
            % (usually the deepest samples).
            
            % Scalar data
            smoothFactor = 0.35; % smoothFactor multiplies number of data points to give loess smoothing parameter
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
                % Variability in model outputs may result in a
                % noisy/non-smooth cost function.
                % Smoothing model outputs before comparing to data may help
                % by improving the signal:noise ratio... fitting to pattern
                % in data rather than each individual data point.
                % Smoothing will also reduce effects of overfitting to data
                % points the model cannot replicate.
                
                smoothF = smoothFactor * n;
                xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothF, 'loess');
                
%                 plot(xx, CDFmod(:,2))
                
                
%                 err2 = (CDFobs(:,1) - CDFmod(:,1)) .^ 2; % squared error
                err2 = (CDFobs(:,1) - xx) .^ 2; % squared error, smoothed
                
                avFun = @mean;
%                 avFun = @geomean;
                
                RMS = (avFun(err2)) .^ 0.5; % average absolute error
                RMS = RMS ./ range(CDFobs(:,1)); % scale to get values more in line with Hellinger distance values... this is ad-hoc method, could be improved...
                
                costComponents.(varLabel) = avFun(RMS); % average over trajectory selections

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
            weight_relVsTot = 3; % weighting factor of relative vs total abundance (relative abundance assumed more reliable)
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cta = costComponents.([varLabel '_autotroph_Tot']);
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    cth = costComponents.([varLabel '_heterotroph_Tot']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
                costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
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

        case 'LeastAbsErr_Hellinger'
            
            % There are problems emerging from using CDF values to
            % calculate the scalar cost values. When modelled values are
            % either much larger or smaller than the data range, the cost
            % becomes insensitive... Probably better to use differences
            % between data and model, rather than converting to CDF
            % values... even though it might be difficult to scale for
            % compatibility with the size-data Hellinger distances.
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                
                % Sort observations
                [yobs_sort, o] = sort(yobs);
                CDFobs = [yobs_sort, ...
                    (1:n)' ./ n]; % Empirical CDF of standardised data
                
%                 figure
%                 plot(CDFobs(:,1), CDFobs(:,2)) % empirical CDF of data
                                
                % Reorder modelled values to match sorted data
                ymod_sort = ymod(o,:);
                CDFmod = [ymod_sort, CDFobs(:,2)]; % align modelled values with the smoothed CDF of the data
                
%                 hold on
%                 
%                 scatter(CDFmod(:,1), CDFmod(:,2))
%                 for ij = 1:n
%                     plot([CDFobs(ij,1), CDFmod(ij,1)], [CDFobs(ij,2), CDFmod(ij,2)], 'Color', [0.85, 0.85, 0.85])
%                 end
                
%                 smoothFactor = 0.35 * n;
%                 xx = smooth(CDFmod(:,2), CDFmod(:,1), smoothFactor, 'loess');                
%                 plot(xx, CDFmod(:,2))

                absErr = abs(CDFobs(:,1) - CDFmod(:,1));
%                 absErr = abs(CDFobs(:,1) - xx);
                
                absErr = absErr ./ range(CDFobs(:,1)); % scale by data range to make cost values more compatible with size data cost

                % Average over data points and trajectory selections
%                 avFun = @mean;
                avFun = @geomean;
%                 avFun = @median;
                avDist = avFun(absErr);  % average over data points
                avFun = @mean;
                costComponents.(varLabel) = avFun(avDist); % average over trajectory selections
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
            
            %             % Within each size data group, weight relative abundance-at-size
            %             % relative to total abundance.
            %             weight_relVsTot = 3; % weighting factor of relative vs total abundance
            %             costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            %             for i = 1:length(waterMasses)
            %                 wm = waterMasses{i};
            %                 if groupedByWaterOrigin
            %                     cta = costComponents.([varLabel '_' wm '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_' wm '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
            %                 else
            %                     cta = costComponents.([varLabel '_autotroph_Tot']);
            %                     cra = costComponents.([varLabel '_autotroph_Rel']);
            %                     cth = costComponents.([varLabel '_heterotroph_Tot']);
            %                     crh = costComponents.([varLabel '_heterotroph_Rel']);
            %                 end
            %                 costSize(1,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cta, cra]);
            %                 costSize(2,i) = mean(2 .* [1, weight_relVsTot] ./ (weight_relVsTot+1) .* [cth, crh]);
            %             end
            
            
            % Omit the total abundance info from the size data --
            % simplified compared to the commented code above...
            costSize = zeros(2,length(waterMasses)); % store weighted costs for all size data groups
            for i = 1:length(waterMasses)
                wm = waterMasses{i};
                if groupedByWaterOrigin
                    cra = costComponents.([varLabel '_' wm '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_' wm '_heterotroph_Rel']);
                else
                    cra = costComponents.([varLabel '_autotroph_Rel']);
                    crh = costComponents.([varLabel '_heterotroph_Rel']);
                end
                costSize(1,i) = cra;
                costSize(2,i) = crh;
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
            %             cost = cost(2); % try fitting only to the size data

            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'meanCDFdist_HellingerFullSpectrum'
            
            % As above, but use model output representing the full size
            % spectra rather than binned (integrated) data.
            % In fact, it is potentially far better than the above because
            % it is now much easier to separately fit all the individual
            % size spectra rather than averaging over sampling events and
            % depths... of course, we could still take averages...
            % I don't know which will be the more useful approach. Fitting
            % to each spectra separately is quite likely to produce poor
            % fits, but the average may still be reasonable. Fitting to the
            % different depth layers could also prove very beneficial by
            % much more strongly informing dependence of size structure
            % upon depth.
            
            % Build the size data cost section to separately fit all size
            % spectra -- keep the returned cost components as autotrophs &
            % heterotrophs in Arctic and/or Atlantic waters => average the
            % costs over sampling events and depths within this function.
            
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                %                 m = size(ymod, 2);
                cdf = (1:n)' ./ n;
                [yobs_sort, o] = sort(yobs); % sort observations
                % put modelled output in same order -- these values will
                % almost certainly NOT be in ascending order...
                ymod_sort = ymod(o,:);
                %                 figure(1)
                %                 plot(yobs_sort, cdf)
                %                 hold on
                %                 plot(ymod_sort, cdf)
                %                 hold off
                % find distances between each modelled data point and the
                % empirical data CDF
                ymodi = interp1(yobs_sort, cdf, ymod_sort, 'linear', 'extrap');
                ymodi = min(1, max(0, ymodi));
                CDFdist = abs(cdf - ymodi);
                avDist = mean(CDFdist); % average over data points
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end
            
            % Vector (size) data
            trophicLevels = unique(Data.size.trophicLevel);
            groupedByWaterOrigin = isfield(Data.size, 'waterMass');
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.waterMass);
            else
                waterMasses = {[]};
            end
            % The VarsSize labels are wrong => redefine
            switch Data.size.obsInCostFunction{1}
                case 'BioVol'; VarsSize = {'BioVolDensity'};
                case 'CellConc'; VarsSize = {'cellDensity'};
            end
            allEvents = unique(Data.size.Event, 'stable');
            nVars = length(VarsSize);
            nWaterMasses = length(waterMasses);
            nTrophicLevels = length(trophicLevels);
            HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
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
                                % Fit only to relative abundance-at-size =>
                                % normalise the data & modelled output
                                yobs = yobs ./ sum(yobs);
                                ymod = ymod ./ sum(ymod);
                                HellingerDistance = (1 - sum((yobs .* ymod) .^ 0.5)) .^ 0.5;
                                HellingerDistances(i,w,it,counter) = HellingerDistance;
                            end
                        end
                    end
                end
            end
            
            % Average Hellinger distances over sample events & depths
            HellingerDistances_ = mean(HellingerDistances, ndims(HellingerDistances));
            % Store in costComponents
            for i = 1:length(VarsSize)
                varLabel = VarsSize{i};
                for it = 1:length(trophicLevels)
                    trophicLevel = trophicLevels{it};
                    if groupedByWaterOrigin
                        for w = 1:length(waterMasses)
                            waterMass = waterMasses{w};
                            label = [varLabel '_' waterMass '_' trophicLevel];
                            costComponents.(label) = HellingerDistances_(i,w,it);
                        end
                    else
                        label = [varLabel '_' trophicLevel];
                        costComponents.(label) = HellingerDistances_(i,1,it);
                    end
                end
            end
            
            % Calculate total cost -- a scalar value
            fields = fieldnames(costComponents);
            costComponents_ = struct2array(costComponents);
            costNutrient = costComponents_(ismember(fields, Vars));
            costSize = costComponents_(~ismember(fields, Vars));
            cost_ = [mean(costNutrient), mean(costSize)];
            cost = mean(cost_);
            
            
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        case 'meanCDFdist_HellingerFullSpectrum_averagedEventsDepths'
            
            % As above, but average over all events and depths before
            % comparing model to data
            
            
            % Scalar data
            for i = 1:length(Vars)
                varLabel = Vars{i};
                yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
                ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
                n = size(ymod, 1);
                %                 m = size(ymod, 2);
                cdf = (1:n)' ./ n;
                [yobs_sort, o] = sort(yobs); % sort observations
                ymod_sort = ymod(o,:);
                % find distances between each modelled data point and the
                % empirical data CDF
                ymodi = interp1(yobs_sort, cdf, ymod_sort, 'linear', 'extrap');
                ymodi = min(1, max(0, ymodi));
                CDFdist = abs(cdf - ymodi);
                avDist = mean(CDFdist); % average over data points
                costComponents.(varLabel) = mean(avDist); % average over trajectory selections
            end
            
            % Vector (size) data
            trophicLevels = unique(Data.size.trophicLevel);
            groupedByWaterOrigin = isfield(Data.size, 'waterMass');
            if groupedByWaterOrigin
                waterMasses = unique(Data.size.waterMass);
            else
                waterMasses = {[]};
            end
            % The VarsSize labels are wrong => redefine
            switch Data.size.obsInCostFunction{1}
                case 'BioVol'; VarsSize = {'BioVolDensity'};
                case 'CellConc'; VarsSize = {'cellDensity'};
            end
            allEvents = unique(Data.size.Event, 'stable');
            nVars = length(VarsSize);
            nWaterMasses = length(waterMasses);
            nTrophicLevels = length(trophicLevels);
            nESD = length(unique(Data.size.ESD));
            
            % Extract modelled output and data into arrays with sample
            % event and depth in the trailing dimension
            yobs = nan(nVars, nWaterMasses, nTrophicLevels, nESD);
            ymod = nan(nVars, nWaterMasses, nTrophicLevels, nESD);
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
                                yobs_ = Data.size.(varLabel)(ind3);
                                ymod_ = modData.size.(varLabel)(ind3);
                                % Fit only to relative abundance-at-size =>
                                % normalise the data & modelled output
                                yobs_ = yobs_ ./ sum(yobs_);
                                ymod_ = ymod_ ./ sum(ymod_);
                                yobs(i,w,it,:,counter) = yobs_;
                                ymod(i,w,it,:,counter) = ymod_;
                            end
                        end
                    end
                end
            end
            % Average over sample events and depths
            yobs = mean(yobs, 5);
            ymod = mean(ymod, 5);
            % Find Hellinger distances betwen model output and data for
            % these averaged spectra, for all variables, water masses and
            % trophic levels.
            HellingerDistances = nan(nVars, nWaterMasses, nTrophicLevels);
            for i = 1:nVars
                for w = 1:nWaterMasses
                    for it = 1:nTrophicLevels
                        yobs_ = yobs(i,w,it,:);
                        ymod_ = ymod(i,w,it,:);
                        HellingerDistance = (1 - sum((yobs_ .* ymod_) .^ 0.5)) .^ 0.5;
                        HellingerDistances(i,w,it) = HellingerDistance;
                    end
                end
            end
            
            % Store in costComponents
            for i = 1:nVars
                varLabel = VarsSize{i};
                for it = 1:nTrophicLevels
                    trophicLevel = trophicLevels{it};
                    if groupedByWaterOrigin
                        for w = 1:nWaterMasses
                            waterMass = waterMasses{w};
                            label = [varLabel '_' waterMass '_' trophicLevel];
                            costComponents.(label) = HellingerDistances(i,w,it);
                        end
                    else
                        label = [varLabel '_' trophicLevel];
                        costComponents.(label) = HellingerDistances(i,1,it);
                    end
                end
            end
            
            % Calculate total cost -- a scalar value
            fields = fieldnames(costComponents);
            costComponents_ = struct2array(costComponents);
            costNutrient = costComponents_(ismember(fields, Vars));
            costSize = costComponents_(~ismember(fields, Vars));
            cost_ = [mean(costNutrient), mean(costSize)];
            cost = mean(cost_);

    end
    
end


%%

if ~exist('cost', 'var')
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if ~exist('costComponents', 'var')
    costComponents = nan;
end
