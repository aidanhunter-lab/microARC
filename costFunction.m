function [cost, costComponents, modData, out, auxVars] = ...
    costFunction(pars, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions, varargin)

% Returns a scalar 'cost' describing model misfit to data given parameters
% 'pars' and, optionally, a struct 'costComponents' containing the cost 
% ascribed to each separate data type, all model outputs and the modelled
% equivalents of the data.

% Cost function type may be given by name-value pair in varargin, otherwise
% set as the default
defaultCostType = 'polyLikelihood';
costTypes = {'LSS', 'polyLikelihood', 'polyLikelihood2'}; % possible cost functions coded below
selectFunction = defaultCostType;
if ~isempty(varargin)
    i = strcmp(varargin, 'selectFunction');
    if any(i)
        sf = varargin{find(i)+1};
        if ismember(sf, costTypes), selectFunction = sf; end
    end
end

%% Run model

% Set parameter values
Params = updateParameters(Params, FixedParams, pars);

% Integrate
[out, auxVars] = integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);

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
    itraj = evTraj(:,i); % index trajectories used for event
    vars = unique(Data.scalar.Variable(Data.scalar.Event == i)); % variables measured during this event
    vars = vars(ismember(vars, Vars));
    Yearday = Data.scalar.Yearday(find(iEvent, 1));
    Depth = Data.scalar.Depth(iEvent);
    Variable = Data.scalar.Variable(iEvent);
    modData.scalar.Yearday(iEvent) = Yearday;
    modData.scalar.Depth(iEvent) = Depth;
    modData.scalar.Variable(iEvent) = Variable;
    for j = 1:length(vars)
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
et = evTraj(:,ev);
nevent = size(et, 2);
% sample times of each event
etime = nan(nevent, 1);
for i = 1:nevent
    etime(i) = unique(Data.scalar.Yearday(Data.scalar.Event == ev(i)));
end
depths_obs = [unique(Data.size.DepthMin) unique(Data.size.DepthMax)]; % sample depth range
depth_ind = -FixedParams.zw(2:end) > depths_obs(1,1) & ...
    -FixedParams.zw(1:end-1) < depths_obs(1,2); % depth layers corresponding to samples

VarsSize = Data.size.obsInCostFunction;
allVarsSize = unique(Data.size.dataBinned.Variable);

for i = 1:length(allVarsSize)
    vs = allVarsSize{i};
%     vs = VarsSize{i};
    ind = strcmp(Data.size.dataBinned.Variable, vs);
    modData.size.Variable(ind,:) = Data.size.dataBinned.Variable(ind);
    for j = 1:nsamples
        itraj = et(j,:); % trajectory associated with each sampling event
        [~, J] = max(sum(out.P(:,depth_ind,FixedParams.PP_Chl_index,etime,itraj))); % depth layer of modelled chl maximum
        J = squeeze(J);
        switch vs
            case 'NConc'
                ymod = squeeze(out.P(:,depth_ind,FixedParams.PP_N_index,etime,itraj));
                ymod_ = nan(FixedParams.nPP_size, nevent);
                for k = 1:nevent
                    ymod_(:,k) = ymod(:,J(k,k),k,k);
                end
                ymod = mean(ymod_, 2); % average size-spectra over sampling events
            case 'CellConc'
                ymod = auxVars.cellDensity(:,depth_ind,etime,itraj);
                ymod_ = nan(FixedParams.nPP_size, nevent);
                for k = 1:nevent
                    ymod_(:,k) = ymod(:,J(k,k),k,k);
                end
                ymod = mean(ymod_, 2); % average size-spectra over sampling events
        end        
        modData.size.Value(ind,j) = ymod;
        modData.size.scaled_Value(ind,j) = Data.size.(['scaleFun_' vs])( ...
            Data.size.dataBinned.scale_mu(ind), ... 
            Data.size.dataBinned.scale_sig(ind), ymod);
    end
end




%% Cost function

switch selectFunction
    case 'LSS'
        % Least sum of squares
        
        for i = 1:length(Vars)
            varLabel = Vars{i};
            ind = strcmp(modData.scalar.Variable, varLabel);
            ymod = modData.scalar.scaled_Value(ind,:);
            yobs = Data.scalar.scaled_Value(ind);
            y = (ymod - yobs) .^ 2;
            L.(varLabel) = sum(y(:)) / numel(y);
        end
        
        ymod = modData.size.scaled_Ntot;
        yobs = Data.size.dataBinned.scaled_Ntot;
        y = (ymod - yobs) .^ 2;
        L.size = sum(y(:)) / size(y, 2);
        
        costComponents = L;
        cost = sum(struct2array(costComponents));
        
    case 'polyLikelihood'        
        % Model misfit to data described using a 'synthetic likelihood',
        % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
        % The standardised data are represented by polynomial coefficients.
        % Equivalent polynomial coefficients representing modelled output are
        % compared to data in a Gaussian likelihood function. Running the model
        % over multiple forcing data trajectories generates the output variability.
        
        % Fit polynomials to model output
        maxDegree = Data.scalar.maxDegree;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            ind = strcmp(modData.scalar.Variable, Vars{i});
            x = Data.scalar.(['polyXvals_' varLabel]);
            y = modData.scalar.scaled_Value(ind,:);
            o = Data.scalar.(['sortOrder_' varLabel]); % sorting order
            ys = y(o,:); % model output in same order as sorted data
            polyCoefs = nan(1 + maxDegree, nsamples);
            for j = 1:nsamples % treat each trajectory as random sample of model output specific to each sampling event
                polyCoefs(:,j) = polyfit(x, ys(:,j), maxDegree);
            end
            modData.scalar.(['polyCoefs_' varLabel]) = polyCoefs;
        end
        
        maxDegree = Data.size.maxDegree;
        x = Data.size.polyXvals;
        y = modData.size.scaled_Ntot;
        o = Data.size.sortOrder;
        ys = y(o,:); % model output in same order as sorted data
        polyCoefs = nan(1 + maxDegree, nsamples);
        for j = 1:nsamples % treat each trajectory as random sample of model output specific to each sampling event
            polyCoefs(:,j) = polyfit(x, ys(:,j), maxDegree);
        end
        modData.size.polyCoefs = polyCoefs;
                
        % Cost function
        % Negative log-likelihood based on polynomial coefficients describing data
        % shape. Distributions of coefficients describing model output are
        % generated by running the model over multiplevforcing data trajectories.
        log2pi = log(2*pi);
        for i = 1:length(Vars)
            varLabel = Vars{i};
            coefs = modData.scalar.(['polyCoefs_' varLabel]);
            ncoefs = size(coefs, 1);
            mu = mean(coefs, 2);
            S = coefs - mu;
            sig = (S * S') / (size(S, 2) - 1);
            sig = 0.5 * (sig + sig'); % guarentees symmetry
            if ~all(eig(sig) > 0) % if not positive definite then coerce sig to SPD
                sig = nearestSPD(sig);
            end
            y = Data.scalar.(['polyCoefs_' Vars{i}])(:) - mu;
            sigChol = chol(sig);             % numerically stable determinant of
            logdetsig = 2*sum(log(diag(sigChol))); % covariance matrix
            L.(Vars{i}) = 0.5 * (y' * (sig \ y) + logdetsig + ncoefs * log2pi);
        end
        
        coefs = modData.size.polyCoefs;
        ncoefs = size(coefs, 1);
        mu = mean(coefs, 2);
        S = coefs - mu;
        sig = (S * S') / (size(S, 2) - 1);
        sig = 0.5 * (sig + sig');
        if ~all(eig(sig) > 0)
            sig = nearestSPD(sig);
        end
        y = Data.size.polyCoefs(:) - mu;
        sigChol = chol(sig);
        logdetsig = 2*sum(log(diag(sigChol)));
        L.size = 0.5 * (y' * (sig \ y) + logdetsig + ncoefs * log2pi);
        
        costComponents = L;
        cost = sum(struct2array(costComponents));
    
    case 'polyLikelihood2' % include some weighting factors for polynomial coefficients
        % Model misfit to data described using a 'synthetic likelihood',
        % as described by Wood (2010), Nature Letters, 466. doi:10.1038/nature09319
        % The standardised data are represented by polynomial coefficients.
        % Equivalent polynomial coefficients representing modelled output are
        % compared to data in a Gaussian likelihood function. Running the model
        % over multiple forcing data trajectories generates the output variability.
        
        % Fit polynomials to model output
        maxDegree = Data.scalar.maxDegree;
        for i = 1:length(Vars)
            varLabel = Vars{i};
            ind = strcmp(modData.scalar.Variable, Vars{i});
            x = Data.scalar.(['polyXvals_' varLabel]);
            y = modData.scalar.scaled_Value(ind,:);
            o = Data.scalar.(['sortOrder_' varLabel]); % sorting order
            ys = y(o,:); % model output in same order as sorted data
            polyCoefs = nan(1 + maxDegree, nsamples);
            for j = 1:nsamples % treat each trajectory as random sample of model output specific to each sampling event
                polyCoefs(:,j) = polyfit(x, ys(:,j), maxDegree);
            end
            modData.scalar.(['polyCoefs_' varLabel]) = polyCoefs;
        end
        
        maxDegree = Data.size.maxDegree;
        for i = 1:length(allVarsSize)
            varLabel = allVarsSize{i};
            ind = strcmp(modData.size.Variable, varLabel);
            x = Data.size.(['polyXvals_' allVarsSize{i}]);
            y = modData.size.scaled_Value(ind,:);
            o = Data.size.(['sortOrder_' varLabel]);
            ys = y(o,:); % model output in same order as sorted data
            polyCoefs = nan(1 + maxDegree, nsamples);
            for j = 1:nsamples % treat each trajectory as random sample of model output specific to each sampling event
                polyCoefs(:,j) = polyfit(x, ys(:,j), maxDegree);
            end
            modData.size.(['polyCoefs_' varLabel]) = polyCoefs;
        end
        
        % Cost function
        % Use uni-variate normal likelihoods separately for each polynomial
        % coefficient, weight the different terms, then combine into single
        % likelihood. This neglects covariance between coefficients...
        log2pi = log(2*pi);
        for i = 1:length(Vars)            
            varLabel = Vars{i};
            coefs = modData.scalar.(['polyCoefs_' varLabel]);
            mu = mean(coefs, 2);
            S = coefs - mu;
            Cov = (S * S') / (size(S, 2) - 1);
            Cov = 0.5 * (Cov + Cov'); % guarentees symmetry
            if ~all(eig(Cov) > 0) % if not positive definite then coerce sig to SPD
                Cov = nearestSPD(Cov);
            end            
            V = diag(Cov); % variance of modelled polynomial coefficients
            y = Data.scalar.(['polyCoefs_' varLabel])(:) - mu;
            weights = Data.scalar.(['polyWeights_' varLabel]);
            weights = weights ./ sum(weights) .* length(weights); % rescale weights about 1
            negLogLik = 0.5 * (log(V) + y .^ 2 ./ V + log2pi);
            weightedNegLogLik = weights(:) .* negLogLik;
            L.(varLabel) = sum(weightedNegLogLik) / length(y);
        end
        
        for i = 1:length(VarsSize)
            varLabel = VarsSize{i};
            coefs = modData.size.(['polyCoefs_' varLabel]);
            mu = mean(coefs, 2);
            S = coefs - mu;
            Cov = (S * S') / (size(S, 2) - 1);
            Cov = 0.5 * (Cov + Cov');
            if ~all(eig(Cov) > 0)
                Cov = nearestSPD(Cov);
            end
            V = diag(Cov);
            y = Data.size.(['polyCoefs_' varLabel])(:) - mu;
            weights = Data.size.(['polyWeights_' varLabel]);
            weights = weights ./ sum(weights) .* length(weights);
            negLogLik = 0.5 * (log(V) + y .^ 2 ./ V + log2pi);
            weightedNegLogLik = weights(:) .* negLogLik;
            L.(varLabel) = sum(weightedNegLogLik) / length(y);
        end
        
        costComponents = L;
        cost = sum(struct2array(costComponents));

        
end
