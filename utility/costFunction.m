function [cost, costComponents, modData, out, auxVars] = costFunction(pars, ... 
    FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, varargin)
% Returns a scalar 'cost' quantifying model misfit to data, given parameters
% 'pars'.
% May optionally return a struct 'costComponents' containing the cost 
% ascribed to each separate data type; all model outputs, 'out' and 
% 'auxVars'; and the modelled equivalents of the data, 'modData'.

extractVarargin(varargin)

returnExtra_ = {'cellDensity', 'biovolume'}; % extra output needed for the cost function
if ~exist('returnExtra', 'var')
    % more outputs may be returned as extra outputs (not recommended while
    % optimising parameters)
    returnExtra = returnExtra_;
else
    returnExtra = unique([eval('returnExtra'), returnExtra_], 'stable');
end
    
% Cost function type should be given by name-value pair in varargin,
% otherwise the default is used
defaultCostType = 'LSS';
if ~exist('selectFunction', 'var')
    selectFunction = defaultCostType;
end

Vars = Data.scalar.obsInCostFunction; % Data types used to fit the model
VarsSize = Data.size.obsInCostFunction;


%% Run model

% Set parameter values
Params = updateParameters(Params, FixedParams, pars);

% % Some initial state variables are specified as functions of parameter
% % values => should be recalculated for each parameter set. If v0 were totally
% % independent of parameter values then it would be absolutely fixed and, to 
% % avoid wasted computation, should be passed as an argument to costFunction.m.
% v0 = initialiseVariables(FixedParams, Params, Forc);

% Integrate
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ...
    odeIntegrator, odeOptions, 'returnExtra', returnExtra);


%% Match model outputs to data

% Extract output from times and depths matching the data, then transform
% the output using the same functions that standardised the data. The
% resulting 'modData' should be comparable to the observations in 'Data'.
modData = matchModOutput2Data(out, auxVars, Data, FixedParams);


%% Cost function

switch selectFunction
    case 'LSS'
        % Least sum of squares
        
        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};            
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);            
            squaredError = (yobs - ymod) .^ 2;
            L.(varLabel) = sum(squaredError(:)) / numel(squaredError);
        end
        
        % Size spectra data
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            % autotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            yobs = Data.size.dataBinned.scaled_Value(ind);
            ymod = modData.size.scaled_Value(ind,:);
            squaredError = (yobs - ymod) .^ 2;
            L.([varLabel '_autotroph']) = sum(squaredError(:)) / numel(squaredError);
            % heterotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
            yobs = Data.size.dataBinned.scaled_Value(ind);
            ymod = modData.size.scaled_Value(ind,:);
            squaredError = (yobs - ymod) .^ 2;
            L.([varLabel '_heterotroph']) = sum(squaredError(:)) / numel(squaredError);
            L.(varLabel) = L.([varLabel '_autotroph']) + L.([varLabel '_heterotroph']);
        end
        
        costComponents = L;
        cost = 0;
        for i = 1:length(Vars)
            cost = cost + costComponents.(Vars{i});
        end
%         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
        for i = 1:length(VarsSize)
            cost = cost + costComponents.(VarsSize{i});
        end
        
        
    case 'RMS'
        % Root mean square
        
        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            absError = sqrt((yobs - ymod) .^ 2);
            L.(varLabel) = sum(absError(:)) / numel(absError);
        end
        
        % Size spectra data
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            % autotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');            
            yobs = Data.size.dataBinned.scaled_Value(ind);
            ymod = modData.size.scaled_Value(ind,:);
            absError = sqrt((yobs - ymod) .^ 2);
            L.([varLabel '_autotroph']) = sum(absError(:)) / numel(absError);
            % heterotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');            
            yobs = Data.size.dataBinned.scaled_Value(ind);
            ymod = modData.size.scaled_Value(ind,:);
            absError = sqrt((yobs - ymod) .^ 2);
            L.([varLabel '_heterotroph']) = sum(absError(:)) / numel(absError);
            L.(varLabel) = L.([varLabel '_autotroph']) + L.([varLabel '_heterotroph']);
        end
        
        costComponents = L;
        cost = 0;
        for i = 1:length(Vars)
            cost = cost + costComponents.(Vars{i});
        end
        %         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
        for i = 1:length(VarsSize)
            cost = cost + costComponents.(VarsSize{i});
        end
        
    case 'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormalDirichlet' % fit every data point...
        % Model misfit to data described using a 'synthetic likelihood',
        % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
        % Standardised scalar data are approximately normally distributed.
        % A variety of trajectory combinations are used to generate model
        % outputs, which are transformed identically to the data. Model
        % outputs define normal distributions used as likelihood terms --
        % what is likelihood of observing data given model expectations?
        % Running the model over multiple forcing data trajectories
        % generates the output variability (process error which is probably
        % underestimated).
        
        % Scalar data
        % Each standardised data point is assigned an independent normal 
        % distribution parameterised using model outputs over multiple 
        % trajectories. The likelihood is the product of probabilities of 
        % observing all data points given these model-estimated 
        % distributions.

        log2pi = log(2*pi);
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            mu = mean(ymod, 2); % expectation of each data point
            sig2 = var(ymod, 0, 2); % variance of model output
            sig2(sig2 == 0) = min(sig2(sig2 > 0)); % include for robustness (we only see zero variability when using single trajectories for any sampling event)
            n = length(yobs); % sample size
%             L = prod(1 ./ ((2*pi*sig2) .^ 0.5) .* exp(-0.5 ./ sig2 .* (yobs - mu) .^ 2)) .^ (1/n);
            negLogLik = 0.5 .* (log2pi + 1/n .* sum(log(sig2) + (yobs - mu) .^ 2 ./ sig2));
            L.(varLabel) = negLogLik;
        end

        % Size data        
        % Dirichlet distribution on relative abundance info in size
        % spectra. Log-normal distribution on total abundance in size data        
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            
            % autotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            yobs = Data.size.dataBinned.Value(ind);
            ymod = modData.size.Value(ind,:);
            % Derive Dirichlet distribution parameters from model output.
            yobsTot = sum(yobs);
            ymodTot = sum(ymod);
            pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
            pmod = ymod ./ ymodTot;
            alpha = fitDirichlet(pmod); % estimate concentration parameter
            % Dirichlet likelihood for simplex
%            L = gamma(alpha0) ./ prod(gamma(alpha)) .* prod(p_obs .^ (alpha-1));
            negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));            
            
            % Lognormal likelihood for total
            yobsTot_log = log(yobsTot);
            ymodTot_log = log(ymodTot);
            mu = mean(ymodTot_log);
            sig2 = var(ymodTot_log);
            negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
            L.([varLabel '_Rel_autotroph']) = negLogLik;
            L.([varLabel '_Tot_autotroph']) = negLogLik2;

            % heterotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
            yobs = Data.size.dataBinned.Value(ind);
            ymod = modData.size.Value(ind,:);
            
            if FixedParams.nZP_size ~= 1
                
                % Derive Dirichlet distribution parameters from model output.
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
                pmod = ymod ./ ymodTot;
                alpha = fitDirichlet(pmod); % estimate concentration parameter
                % Dirichlet likelihood for simplex
                %            L = gamma(alpha0) ./ prod(gamma(alpha)) .* prod(p_obs .^ (alpha-1));
                negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));
                
                % Lognormal likelihood for total
                yobsTot_log = log(yobsTot);
                ymodTot_log = log(ymodTot);
                mu = mean(ymodTot_log);
                sig2 = var(ymodTot_log);
                negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
                L.([varLabel '_Rel_heterotroph']) = negLogLik;
                L.([varLabel '_Tot_heterotroph']) = negLogLik2;
                
                L.(varLabel) = L.([varLabel '_Rel_autotroph']) + L.([varLabel '_Tot_autotroph']) + ... 
                    L.([varLabel '_Rel_heterotroph']) + L.([varLabel '_Tot_heterotroph']);

            else
                
                yobsTot = sum(yobs);
                ymodTot = ymod(1,:);

                % Lognormal likelihood for total
                yobsTot_log = log(yobsTot);
                ymodTot_log = log(ymodTot);
                mu = mean(ymodTot_log);
                sig2 = var(ymodTot_log);
                negLogLik = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
                L.([varLabel '_Tot_heterotroph']) = negLogLik;
                
                L.(varLabel) = L.([varLabel '_Rel_autotroph']) + L.([varLabel '_Tot_autotroph']) + ...
                    L.([varLabel '_Tot_heterotroph']);

            end
            
        end
        
        costComponents = L;
        cost = 0;
        for i = 1:length(Vars)
            cost = cost + costComponents.(Vars{i});
        end
%         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
        for i = 1:length(VarsSize)
            cost = cost + costComponents.(VarsSize{i});
        end
%         cost = sum(struct2array(costComponents));
    
    
    case 'syntheticLikelihood_ScalarNormal_SizeSpectraLogNormal_logisticNormal'
        
        % Scalar data
        
        log2pi = log(2*pi);
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1); % sample size
            m = size(ymod, 2); % number of replicates
            mu = mean(ymod); % modelled expectations and 
            sig = std(ymod); % standard deviations for selected trajectories
            sig2 = sig .^ 2;
            negLogLik = 0.5 .* (log2pi + log(sig2) + (yobs - mu) .^ 2 ./ sig2);
            negLogLik = sum(negLogLik(:)) / n / m;
            L.(varLabel) = negLogLik;
        end
        
        % Size data
        
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.size.dataBinned.Variable, varLabel);
            
            % autotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            
            yobs = Data.size.dataBinned.Value(ind);
            ymod = modData.size.Value_allEvents(ind,:,:); % modelled values for each separate sampling event
            ymodMean = modData.size.Value(ind,:); % modelled values averaged over sampling events
            % dimension of ymod is [size, sampling event, replicate (trajectory choice), sample number (year or cruise)]
            nsize = size(ymod, 1); % number of sizes
%             ne = size(ymod, 2); % number of sampling events
            m = size(ymod, 3); % number of replicates (trajectory choices)
%             n = size(ymod, 4); % number of samples/years/cruises etc
            
            % decompose spectra into total and relative abundance (a single number and a simplex)
            yobsTot = sum(yobs); % totals
            ymodTot = sum(ymod);
            ymodMeanTot = sum(ymodMean);
            pobs = yobs ./ yobsTot; % simplices
            pmod = ymod ./ ymodTot;
            pmodMean = ymodMean ./ ymodMeanTot;
            
            % Use a logistic-normal distribution to calculate likelihood of
            % observed relative abundance (simplex) given the modelled
            % equivalents.
            % Methods detailed in: Francis, R.I.C.C. (2014) Fish.Res.151:70-84. doi:10.1016/j.fishres.2013.12.015
            
            % Assuming observed simplex, pobs, has a logistic-normal
            % distribution with expectation, pmod => X is multi-variate normal
            % if pobs = exp(X) / sum(exp(X)), where X has expectation, log(pmod),
            % and covariance, C.
            % The pobs-X transform is not one-to-one, so transform pobs
            % by reducing the size dimension by one.
            Yobs = log(pobs(1:end-1) ./ pobs(end));
            % Now assume that Yobs is multi-variate normal with
            % expectation, mu, and covariance, V.
%             mu = log(pmod(1:end-1,:,:) ./ pmod(end,:,:)); % expectations for [nsize-1] multi-variate normal Yobs (all sampling events)
            Mu = log(pmodMean(1:end-1,:) ./ pmodMean(end,:)); % expectations for [nsize-1] multi-variate normal Yobs (averaged sampling events)

            log_pmod = log(pmod); % expectations (for each separate sampling event) of [nsize] multivariate-normal distribution for X
            
            K = [eye(nsize-1) -ones(nsize-1,1)];
            
            % Use modelled output for each separate sampling event to
            % estimate (co)variance parameter for the event-averaged data
            
            C = nan(nsize, nsize, m); % Covaraince of multi-variate normal distribution
            V = nan(nsize-1, nsize-1, m); % Transformed covariance -- useful form for likelihood
            
            covarianceTypes = {'var', 'CVtridiag', 'CVfull'}; % variance, tri-diagonal covariance matrix, full covariance matrix
            covarianceType = covarianceTypes{1};
            
            switch covarianceType
                case 'var'
                    for j = 1:m
                        C(:,:,j) = diag(diag(cov(log_pmod(:,:,j)')));
                        V(:,:,j) = K * C(:,:,j) * K';
                    end
                case 'CVfull'
                    for j = 1:m
                        C(:,:,j) = cov(log_pmod(:,:,j)');
                        V(:,:,j) = K * C(:,:,j) * K';
                    end
                case 'CVtridiag'
                    for j = 1:m
                        C(:,:,j) = cov(log_pmod(:,:,j)');
                        C(:,:,j) = diag(diag(C(:,:,j))) + ... 
                            diag(diag(C(:,:,j), -1), -1) + ... 
                            diag(diag(C(:,:,j), 1), 1);
                        V(:,:,j) = K * C(:,:,j) * K';
                    end
            end
            
            w = Yobs - Mu; % error vectors
            
            for j = 1:m
                negLogLik(j) = log(det(V(:,:,j))) + (w(:,j)' / V(:,:,j) * w(:,j));
            end
%             negLogLik = 0.5 .* (negLogLik + (nsize - 1) .* log2pi) + sum(log(pobs));
            negLogLik = 0.5 .* (negLogLik + (nsize - 1) .* log2pi);
            negLogLik_rel = sum(negLogLik) / m;
%             negLogLik_rel = sum(negLogLik) / m / (nsize-1);
            
            % Lognormal likelihood for total abundance
            yobsTot_log = log(yobsTot);
            ymodTot_log = log(ymodTot);
            mu = squeeze(mean(ymodTot_log));
            sig = squeeze(std(ymodTot_log));
            sig2 = sig .^ 2;
            
%             Lik = prod(1 ./ ((2*pi*sig2) .^ 0.5) .* exp(-0.5 .* (yobsTot_log - mu) .^ 2 ./ sig2)) .^ (1/m)            
            negLogLik_tot = 0.5 .* (m .* log2pi + sum(log(sig2) + (yobsTot_log - mu) .^ 2 ./ sig2)) ./ m;
            
            L.([varLabel '_Rel_autotroph']) = negLogLik_rel;
            L.([varLabel '_Tot_autotroph']) = negLogLik_tot;
            L.([varLabel '_autotroph']) = negLogLik_rel + negLogLik_tot;

            
            % heterotrophs
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'heterotroph');
            
            yobs = Data.size.dataBinned.Value(ind);
            ymod = modData.size.Value_allEvents(ind,:,:); % modelled values for each separate sampling event
%             ymodMean = modData.size.Value(ind,:); % modelled values averaged over sampling events
            
            singleSizeClass = FixedParams.nZP_size == 1;
            
            if singleSizeClass
                ymod = ymod(1,:,:);
%                 ymodMean = ymodMean(1,:);
                yobs = sum(yobs); % sum observations over size classes
                
                % dimension of ymod is [size, sampling event, replicate (trajectory choice), sample number (year or cruise)]
                m = size(ymod, 3); % number of replicates (trajectory choices)
                
                % Lognormal likelihood for total abundance
                yobsTot_log = log(yobs);
                ymodTot_log = log(ymod);
                mu = squeeze(mean(ymodTot_log));
                sig2 = squeeze(var(ymodTot_log));
                
                negLogLik_tot = 0.5 .* (m .* log2pi + sum(log(sig2) + (yobsTot_log - mu) .^ 2 ./ sig2)) ./ m;

                L.([varLabel '_heterotroph']) = negLogLik_tot;
                
            else
                
                warning('still need to set up the multiple size class model!')
                
                
            end
            
            L.(varLabel) = L.([varLabel '_autotroph']) + L.([varLabel '_heterotroph']);
            
        end
        
        costComponents = L;
        cost = 0;
        for i = 1:length(Vars)
            cost = cost + costComponents.(Vars{i});
        end
%         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
        for i = 1:length(VarsSize)
            cost = cost + costComponents.(VarsSize{i});
        end
%         cost = sum(struct2array(costComponents));



    case 'syntheticLikelihood_ScalarNormalShape_SizeSpectraLogNormalDirichlet'
        % Model misfit to data described using a 'synthetic likelihood',
        % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
        % Standardised scalar data are approximately normally distributed.
        % A variety of trajectory combinations are used to generate model
        % outputs, which are transformed identically to the data. Model
        % outputs define normal distributions used as likelihood terms --
        % what is likelihood of observing data given model expectations?
        % Running the model over multiple forcing data trajectories
        % generates the output variability (process error which is probably
        % underestimated).
        
        % Scalar data
        % Each standardised data point is assigned an independent normal 
        % distribution parameterised using model outputs over multiple 
        % trajectories. The likelihood is the product of probabilities of 
        % observing all data points given these model-estimated 
        % distributions.
        log2pi = log(2*pi);
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
%             n = length(yobs);
            n = size(ymod, 2);
            mu_obs = mean(yobs);
            v_obs = var(yobs);
            Mu = mean(ymod);
            V = var(ymod);
            % likelihood of mu_obs
            mu = mean(Mu);
            v = var(Mu);
            negLogLik_mu = log2pi + log(v) + (mu_obs - mu) .^ 2 ./ v;
            % likelihood of v_obs
            mu = mean(V);
            v = var(V);
            negLogLik_v = log2pi + log(v) + (v_obs - mu) .^ 2 ./ v;
            % combine
            negLogLik = 0.5 * n * (negLogLik_mu + negLogLik_v);
            L.(varLabel) = negLogLik;
        end
        

        % Size data        
        % Dirichlet distribution on relative abundance info in size
        % spectra. Log-normal distribution on total abundance in size data        
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            yobs = Data.size.dataBinned.Value(strcmp(Data.size.dataBinned.Variable, varLabel));
            ymod = modData.size.Value(strcmp(modData.size.Variable, varLabel),:);
            % Derive Dirichlet distribution parameters from model output.
            yobsTot = sum(yobs);
            ymodTot = sum(ymod);
            pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
            pmod = ymod ./ ymodTot;
            alpha = fitDirichlet(pmod); % estimate concentration parameter
            % Dirichlet likelihood for simplex
%            L = gamma(alpha0) ./ prod(gamma(alpha)) .* prod(p_obs .^ (alpha-1));
            negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));            
            
            % Lognormal likelihood for total
            yobsTot_log = log(yobsTot);
            ymodTot_log = log(ymodTot);
            mu = mean(ymodTot_log);
            sig = std(ymodTot_log);
            sig2 = sig .^ 2;
            negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
            L.([varLabel '_Rel']) = negLogLik;
            L.([varLabel '_Tot']) = negLogLik2;
            L.(varLabel) = negLogLik + negLogLik2;
        end
        
        costComponents = L;
        cost = 0;
        for i = 1:length(Vars)
            cost = cost + costComponents.(Vars{i});
        end
%         cost = cost / length(Vars); % group all scalar variables together, with weighting equal to the size spectra
        for i = 1:length(VarsSize)
            cost = cost + costComponents.(VarsSize{i});
        end
%         cost = sum(struct2array(costComponents));


    case 'N_LN-Dir_groupWaterOrigin'
        % Synthetic likelihoods where variability parameters are estimated
        % using variability in model outputs from the various trajectories.
        % Normal distributions for scalar data points.
        % Size spectra vectors decomposed into totals and simplexes --
        % lognormal distributions for the totals, Dirichlet distributions
        % for the simplexes.
        % Size spectra data for autotrophs and heterotrophs, and for water 
        % of Arctic and of Atlantic origin -- 4 separate data vectors.
        
        % Scalar data
        log2pi = log(2*pi);
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            mu = mean(ymod, 2); % expectation of each data point
            sig2 = var(ymod, 0, 2); % variance of model output
            sig2(sig2 == 0) = min(sig2(sig2 > 0)); % include for robustness (we only see zero variability when using single trajectories for any sampling event)
            n = length(yobs); % sample size
            %             L = prod(1 ./ ((2*pi*sig2) .^ 0.5) .* exp(-0.5 ./ sig2 .* (yobs - mu) .^ 2)) .^ (1/n);
            negLogLik = 0.5 .* (log2pi + 1/n .* sum(log(sig2) + (yobs - mu) .^ 2 ./ sig2));
            L.(varLabel) = negLogLik;
        end
        
        % Vector (size) data
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.sizeFull.dataBinned.groupedByOrigin.Variable, varLabel);
            waterMasses = unique(Data.sizeFull.dataBinned.groupedByOrigin.waterMass);
            
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, wm);
                
                % autotrophs
                ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'autotroph');
                yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
                ymod = modData.sizeFull.(['Value_' wm])(ind,:);
                % Derive Dirichlet distribution parameters from model output.
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
                pmod = ymod ./ ymodTot;
                alpha = fitDirichlet(pmod); % estimate concentration parameter
                % Dirichlet likelihood for simplex
                negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));

                % Lognormal likelihood for total
                yobsTot_log = log(yobsTot);
                ymodTot_log = log(ymodTot);
                mu = mean(ymodTot_log);
                sig2 = var(ymodTot_log);
                negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
                
                L.([varLabel '_' wm '_autotroph_Rel']) = negLogLik;
                L.([varLabel '_' wm '_autotroph_Tot']) = negLogLik2;
                
                % heterotrophs
                ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'heterotroph');
                yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
                ymod = modData.sizeFull.(['Value_' wm])(ind,:);
                
                % Derive Dirichlet distribution parameters from model output.
                yobsTot = sum(yobs);
                ymodTot = sum(ymod);
                pobs = yobs ./ yobsTot; % observed relative abundance -- simplex
                pmod = ymod ./ ymodTot;
                alpha = fitDirichlet(pmod); % estimate concentration parameter
                % Dirichlet likelihood for simplex
                negLogLik = sum(gammaln(alpha)) - gammaln(sum(alpha)) - sum((alpha - 1) .* log(pobs));
                
                % Lognormal likelihood for total
                yobsTot_log = log(yobsTot);
                ymodTot_log = log(ymodTot);
                mu = mean(ymodTot_log);
                sig2 = var(ymodTot_log);
                negLogLik2 = 0.5 .* (log2pi + sum(log(sig2) + 1 ./ sig2 .* (yobsTot_log - mu) .^ 2));
                
                L.([varLabel '_' wm '_heterotroph_Rel']) = negLogLik;
                L.([varLabel '_' wm '_heterotroph_Tot']) = negLogLik2;
            end
        end
        
        costComponents = L;
        
        % Apply any weightings to data types here, before summing
        % costComponents to find total cost...
        cost = sum(struct2array(costComponents));

    
    case 'Hellinger_groupWaterOrigin'
        % Calulate Hellinger distances between observations and their 
        % equivalent model outputs.
        % For consistency, use this metric for all data types.
        
        % second attempt
        
        % Scalar data
        for i = 1:length(Vars)
            varLabel = Vars{i};
            yobs = Data.scalar.scaled_Value(strcmp(Data.scalar.Variable, varLabel));
            ymod = modData.scalar.scaled_Value(strcmp(modData.scalar.Variable, varLabel),:);
            n = size(ymod, 1);
            m = size(ymod, 2);
            % Try calculating cdfs and pdfs directly without assuming
            % normality of model output distributions.
            % To compare observed and modelled distributions, the pdfs will
            % need to be evaluated across identical domains.
            cdf = (1:n)' ./ n;
            yobsc = sort(yobs);                        
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
            % define regular grid across measurement space
            vrange = [min([yobsc(:); ymodc(:)]), max([yobsc(:); ymodc(:)])];
            grid = nan(n, 1);
            grid(2:end) = linspace(vrange(1), vrange(2), n-1)';
            sp = diff(grid(2:3));
            grid(1) = grid(2) - sp;
            % interpolate cdfs and derive pdfs
            cdfobsi = interp1(yobsc, cdfobs, grid, 'linear', 'extrap');
            cdfobsi(cdfobsi < min(cdfobs)) = 0;
            cdfobsi(cdfobsi > 1) = 1;
            pdfobsi = diff(cdfobsi);            
            cdfmodi = nan(n, m);
            pdfmodi = nan(n-1, m);
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
        % If the Hellinger distance metric is used for the size data, and
        % if it requires probabilities as input, then I don't know how to
        % include magnitudes in the cost function. Is this Hellinger
        % distance metric only appropriate for comparing the relative
        % abundance-at-size?
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            ind0 = strcmp(Data.sizeFull.dataBinned.groupedByOrigin.Variable, varLabel);
            waterMasses = unique(Data.sizeFull.dataBinned.groupedByOrigin.waterMass);
            
            for w = 1:length(waterMasses)
                wm = waterMasses{w};
                ind1 = ind0 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.waterMass, wm);
                
                % autotrophs
                ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'autotroph');
                yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
                ymod = modData.sizeFull.(['Value_' wm])(ind,:);
                n = size(ymod, 1);
                m = size(ymod, 2);
                cdfobs = cumsum(yobs) ./ sum(yobs);
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ sum(ymod);
                pdfmod = diff([zeros(1, m); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.([varLabel '_' wm '_autotroph']) = mean(hellingerDistance); % average over trajectory selections

                % heterotrophs
                ind = ind1 & strcmp(Data.sizeFull.dataBinned.groupedByOrigin.trophicLevel, 'heterotroph');
                yobs = Data.sizeFull.dataBinned.groupedByOrigin.Value(ind);
                ymod = modData.sizeFull.(['Value_' wm])(ind,:);                
                n = size(ymod, 1);
                m = size(ymod, 2);
                cdfobs = cumsum(yobs) ./ sum(yobs);
                pdfobs = diff([0; cdfobs]);
                cdfmod = cumsum(ymod) ./ sum(ymod);
                pdfmod = diff([zeros(1, m); cdfmod]);
                hellingerDistance = (1 - sum((pdfobs .* pdfmod) .^ 0.5)) .^ 0.5;
                costComponents.([varLabel '_' wm '_heterotroph']) = mean(hellingerDistance);

            end
        end
        
        % Take averages across data-types to assign equal weightings to the
        % scalar data and size data.
        cost = zeros(1,2);
        fields = fieldnames(costComponents);
        for i = 1:length(fields)
            if ismember(fields{i}, Vars)
                % scalar data
               cost(1) = cost(1) + costComponents.(fields{i});
            else
                % size data
                cost(2) = cost(2) + costComponents.(fields{i});
            end
        end
        cost(1) = cost(1) ./ length(Vars);
        cost(2) = cost(2) ./ (length(fields) - length(Vars));
        cost = mean(cost);
        
%         cost = sum(struct2array(costComponents));

%         disp(costComponents)
%         disp(cost)
        
end


%%

if ~exist('cost', 'var')
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if ~exist('costComponents', 'var')
    costComponents = nan;
end
if ~exist('modData', 'var')
    modData = nan;
end
if ~exist('out', 'var')
    out = nan;
end
if ~exist('auxVars', 'var')
    auxVars = nan;
end

