function [cost, costComponents, modData, out, auxVars] = ...
    costFunction(pars, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, varargin)

% Returns a scalar 'cost' describing model misfit to data given parameters
% 'pars' and, optionally, a struct 'costComponents' containing the cost 
% ascribed to each separate data type, all model outputs and the modelled
% equivalents of the data.

% disp(pars') % useful for debugging
% Cost function type may be given by name-value pair in varargin, otherwise
% set as the default
defaultCostType = 'LSS';
selectFunction = defaultCostType;
returnExtra = {'cellDensity', 'biovolume'}; % extra output needed for the cost function
if ~isempty(varargin)
    i = strcmp(varargin, 'selectFunction');
    if any(i)
        selectFunction = varargin{find(i)+1};
    end
    i = strcmp(varargin, 'returnExtra');
    if any(i)
        returnExtra = varargin{find(i)+1};
    end
end

% % Cost function type may be given by name-value pair in varargin, otherwise
% % set as the default
% defaultCostType = 'LSS';
% costTypes = {'LSS', 'polyLikelihood', 'polyLikelihood2', ... 
%     'syntheticLikelihood_normal_Dirichlet', 'syntheticLikelihood_normalShape_Dirichlet', ...
%     'syntheticLikelihood_normal_multinomialDirichlet'}; % possible cost functions coded below
% selectFunction = defaultCostType;
% if ~isempty(varargin)
%     i = strcmp(varargin, 'selectFunction');
%     if any(i)
%         sf = varargin{find(i)+1};
%         if ismember(sf, costTypes), selectFunction = sf; end
%     end
% end

%% Run model

% Set parameter values
Params = updateParameters(Params, FixedParams, pars);

% Set initial state variable values -- some of which are selected using
% parameter values


% Integrate
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, ...
    odeIntegrator, odeOptions, 'returnExtra', returnExtra);


%% Match model outputs to data

% Extract output from times and depths matching the data, then transform
% the output using the same functions that standardised the data.

% Scalar data
Vars = Data.scalar.obsInCostFunction;
% Vars = unique(Data.scalar.Variable);

evTraj = Data.scalar.evTraj;
nsamples = size(evTraj, 1); % number of trajectories used per sampling event

nEvent = Data.scalar.nEvents;
depths_mod = abs(FixedParams.z);

modData.scalar.Yearday = nan(Data.scalar.nSamples,1);
modData.scalar.Depth = modData.scalar.Yearday;
modData.scalar.Variable = cell(Data.scalar.nSamples,1);
modData.scalar.Value = nan(Data.scalar.nSamples, nsamples);
modData.scalar.scaled_Value = modData.scalar.Value;

% Standardise model output with respect to depth and event using linear mixed models
for i = 1:nEvent
    iEvent = Data.scalar.Event == i; % index event
    itraj = evTraj(:,i); % trajectories used for event i
    vars = unique(Data.scalar.Variable(iEvent)); % variables measured in event i
    vars = vars(ismember(vars, Vars));
    iEvent = ismember(Data.scalar.Variable, vars) & iEvent; % omit unmodelled variables from the event index
    Yearday = Data.scalar.Yearday(find(iEvent, 1));
    Depth = Data.scalar.Depth(iEvent);    
    Variable = Data.scalar.Variable(iEvent);
    modData.scalar.Yearday(iEvent) = Yearday;
    modData.scalar.Depth(iEvent) = Depth;
    modData.scalar.Variable(iEvent) = Variable;
    for j = 1:length(vars)
        % loop through all data types sampled during event i
        jvar = vars{j};
        ind = iEvent & strcmp(Data.scalar.Variable, jvar);
        depth = modData.scalar.Depth(ind);
        switch jvar
            % modelled values [depth,trajectory]
            case 'N'
                ymod = squeeze(out.N(:,:,Yearday,itraj));
            case 'PON'
                ymod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_N_index,Yearday,itraj));
            case 'POC'
                ymod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_C_index,Yearday,itraj));
            case 'chl_a'
                ymod = squeeze(sum(out.P(:,:,FixedParams.PP_Chl_index,Yearday,itraj)));
        end
        ymod = interp1(depths_mod, ymod, depth); % interpolate modelled output to match observation depths
        ymod_scaled = Data.scalar.(['scaleFun_' jvar])(Data.scalar.scale_mu(ind), ...
            Data.scalar.scale_sig(ind), ymod); % scale model output using same functions that scaled the data
        modData.scalar.Value(ind,:) = ymod;
        modData.scalar.scaled_Value(ind,:) = ymod_scaled;
    end
end


% Size spectra

ind = ismember(Data.scalar.Year, unique(Data.size.Year)); % index relavent sampling events
ev = unique(Data.scalar.Event(ind));
et = evTraj(:,ev); % sampling events and associated trajectories
nevent = size(et, 2);
% sample times of each event
etime = nan(nevent, 1);
for i = 1:nevent
    ue = unique(Data.scalar.Yearday(Data.scalar.Event == ev(i)));
    etime(i) = ue(1);
%     etime(i) = unique(Data.scalar.Yearday(Data.scalar.Event == ev(i)));
end
depths_obs = [min(Data.size.DepthMin) max(Data.size.DepthMax)]; % sample depth range
depth_ind = -FixedParams.zw(2:end) > depths_obs(1,1) & ...
    -FixedParams.zw(1:end-1) < depths_obs(1,2); % depth layers corresponding to samples

VarsSize = Data.size.obsInCostFunction;
allVarsSize = unique(Data.size.dataBinned.Variable);

n = size(Data.size.dataBinned.Year, 1);
modData.size.Variable = cell(n,1);
modData.size.trophicLevel = cell(n,1);
modData.size.Value = nan(n,nsamples);
modData.size.Value_allEvents = nan(n,nevent,nsamples);
modData.size.scaled_Value = nan(n,1,nsamples);
modData.size.scaled_Value_allEvents = nan(n,nevent,nsamples);


for i = 1:length(allVarsSize)
    vs = allVarsSize{i};
    ind0 = strcmp(Data.size.dataBinned.Variable, vs);
    modData.size.Variable(ind0,:) = Data.size.dataBinned.Variable(ind0);
    modData.size.trophicLevel(ind0,:) = Data.size.dataBinned.trophicLevel(ind0);
    
    ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
    
    for j = 1:nsamples
        itraj = et(j,:); % trajectories associated with sampling event j
        [~, J] = max(sum(out.P(:,depth_ind,FixedParams.PP_Chl_index,etime,itraj))); % use modelled values from chl-max depth layer to compare to data
        J = squeeze(J);
        switch vs
            case 'NConc'
                ymod = squeeze(out.P(1:1:FixedParams.nPP_size,depth_ind,FixedParams.PP_N_index,etime,itraj));
                ymod_ = nan(FixedParams.nPP_size, nevent);
                for k = 1:nevent
                    ymod_(:,k) = ymod(:,J(k,k),k,k);
                end
                ymod = ymod_;
                ymodMean = mean(ymod, 2); % average size-spectra over sampling events
            case 'CellConc'
                ymod = auxVars.cellDensity(1:1:FixedParams.nPP_size,depth_ind,etime,itraj);
                ymod_ = nan(FixedParams.nPP_size, nevent);
                for k = 1:nevent
                    ymod_(:,k) = ymod(:,J(k,k),k,k);
                end
                ymod = ymod_;
                ymodMean = mean(ymod, 2); % average size-spectra over sampling events
            case 'BioVol'
                ymod = auxVars.biovolume(1:FixedParams.nPP_size,depth_ind,etime,itraj);
                ymod_ = nan(FixedParams.nPP_size, nevent);
                for k = 1:nevent
                    ymod_(:,k) = ymod(:,J(k,k),k,k);
                end
                ymod = ymod_;
                ymodMean = mean(ymod, 2); % average size-spectra over sampling events
        end
        modData.size.Value(ind,j) = ymodMean;
        modData.size.Value_allEvents(ind,:,j) = ymod;

        modData.size.scaled_Value(ind,:,j) = Data.size.(['scaleFun_' vs])( ...
            Data.size.dataBinned.scale_mu(ind), ...
            Data.size.dataBinned.scale_sig(ind), ymodMean);
        modData.size.scaled_Value_allEvents(ind,:,j) = Data.size.(['scaleFun_' vs])( ...
            Data.size.dataBinned.scale_mu(ind), ... 
            Data.size.dataBinned.scale_sig(ind), ymod);
    end
end



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
            yobs = Data.size.dataBinned.scaled_Value(strcmp(Data.size.dataBinned.Variable, varLabel));
            ymod = modData.size.scaled_Value(strcmp(modData.size.Variable, varLabel),:);
            squaredError = (yobs - ymod) .^ 2;
            L.(varLabel) = sum(squaredError(:)) / numel(squaredError);
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
            yobs = Data.size.dataBinned.scaled_Value(strcmp(Data.size.dataBinned.Variable, varLabel));
            ymod = modData.size.scaled_Value(strcmp(modData.size.Variable, varLabel),:);
            absError = sqrt((yobs - ymod) .^ 2);
            L.(varLabel) = sum(absError(:)) / numel(absError);
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
            sig = std(ymod, 0, 2); % standard deviation of model output
            sig(sig == 0) = min(sig(sig > 0)); % include for robustness (we only see zero variability when using single trajectories for any sampling event)
            sig2 = sig .^ 2;
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
            ind = ind0 & strcmp(Data.size.dataBinned.trophicLevel, 'autotroph');
            yobs = Data.size.dataBinned.Value(ind);
            ymod = modData.size.Value_allEvents(ind,:,:); % modelled values for each separate sampling event
            ymodMean = modData.size.Value(ind,:); % modelled values averaged over sampling events
            % dimension of ymod is [size, sampling event, replicate (trajectory choice), sample number (year or cruise)]
            nsize = size(ymod, 1); % number of sizes
            ne = size(ymod, 2); % number of sampling events
            m = size(ymod, 3); % number of replicates (trajectory choices)
            n = size(ymod, 4); % number of samples/years/cruises etc
            
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
            mu = log(pmod(1:end-1,:,:) ./ pmod(end,:,:)); % expectations for [nsize-1] multi-variate normal Yobs (all sampling events)
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
            
            L.([varLabel '_Rel']) = negLogLik_rel;
            L.([varLabel '_Tot']) = negLogLik_tot;
            L.(varLabel) = negLogLik_rel + negLogLik_tot;

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


end


%%

if exist('cost', 'var') == 0
    cost = nan;
    warning('Could not evaluate cost... Check that the name-value pair (selectFunction,costFunctionType) is properly specified and corresponds to a viable option within costFunction.m')
end
if exist('costComponents', 'var') == 0
    costComponents = nan;
end
if exist('modData', 'var') == 0
    modData = nan;
end
if exist('out', 'var') == 0
    out = nan;
end
if exist('auxVars', 'var') == 0
    auxVars = nan;
end
