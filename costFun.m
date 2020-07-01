function cost = costFun(x, FixedParams, Params, Forc, Data, v0, ode45options)

% Returns a scalar describing model misfit to data given parameter
% updates, x

nt = FixedParams.nt;
nz = FixedParams.nz;
nPP = FixedParams.nPP;
nOM = FixedParams.nOM;
nTraj = Forc.nTraj;
nEvent = Data.nEvents;
depths_mod = abs(FixedParams.z);

% Set parameter values
Params = updateParameters(Params, FixedParams, x);

% Integrate
[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, ode45options);

% Extract solutions
out.N = reshape(OUT(FixedParams.IN_index,:,:), [1 nz nt nTraj]);
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP nz nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM nz nt nTraj]);

if sum(nExtra) > 0
    for k = 1:nExtra(1)
        auxVars.(namesExtra{k}) = squeeze(AUXVARS(:,k,:,:));
    end
    for k = 1:nExtra(2)
        auxVars.(namesExtra{k+nExtra(1)}) = squeeze(AUXVARS_2d(:,:,k,:,:));
    end
end


% Compare model to observations
cost_N = nan(1,nEvent);
cost_PON = nan(1,nEvent);
for i = 1:nEvent
    iObs = Data.Event == i; % index data
    time = Data.Yearday(find(iObs,1)); % sample time    
    ti = Data.EventTraj(i,:); % index trajectories
    
    vars = unique(Data.Variable(Data.Event == i));
    
    % inorganic nitrogen
    if any(strcmp(vars,'N'))
        ind = iObs & strcmp('N',Data.Variable);
        y_obs = Data.scaled_Value(ind); % measured values
        depths_obs = Data.Depth(ind); % measurement depths
        y_mod = squeeze(out.N(:,:,time,ti)); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths        
        y_mod = Data.scaleFun_N(Data.scale_mu(ind), Data.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_N(i) = sum(sqErr(:)) / numel(y_mod);
    end
    
    % PON
    if any(strcmp(vars,'PON'))
        ind = iObs & strcmp('PON',Data.Variable);
        y_obs = Data.scaled_Value(ind); % measured values
        depths_obs = Data.Depth(ind); % measurement depths
        y_mod = squeeze(out.OM(FixedParams.POM_index,:,time,ti)); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scaleFun_PON(Data.scale_mu(ind), Data.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_PON(i) = sum(sqErr(:)) / numel(y_mod);
    end    
end

cost = nanmean(cost_N) + nanmean(cost_PON);


