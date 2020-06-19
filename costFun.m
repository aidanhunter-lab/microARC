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

% I think i should scale the misfits by standard deviations measured across
% events

cost_N = nan(1,nEvent);
cost_PON = nan(1,nEvent);
for i = 1:nEvent
    iObs = Data.Event == i; % index data
    time = Data.Yearday(find(iObs,1)); % sample time    
    ti = Data.EventTraj(i,:); % index trajectories
    % inorganic nitrogen    
    N = squeeze(out.N(:,:,time,ti)); % modelled values
    iN = iObs & strcmp('N',Data.Variable);
    N_obs = Data.Value(iN); % measured values
    N_depth = Data.Depth(iN);
    N = interp1(depths_mod, N, N_depth); % interpolate modelled output to match observation depths
    
    cost_N(i) = sum(sum((log(N ./ N_obs)).^2, 2) ./ size(N,2)) ./ size(N,1);

    % PON
    PON = squeeze(out.OM(FixedParams.POM_index,:,time,ti));
    iPON = iObs & strcmp('PON',Data.Variable);
    PON_obs = Data.Value(iPON);
    PON_depth = Data.Depth(iPON);
    keep = PON_depth < max(depths_mod);
    PON_obs = PON_obs(keep);
    PON_depth = PON_depth(keep);
    PON = interp1(depths_mod, PON, PON_depth);

    cost_PON(i) = sum(sum((log(PON ./ PON_obs)).^2, 2) ./ size(PON,2)) ./ size(PON,1);
        
end

cost_N = cost_N(~isnan(cost_N));
cost_PON = cost_PON(~isnan(cost_PON));

cost = mean(cost_N) + mean(cost_PON);





