function [cost, costComponents] = costFun(x, FixedParams, Params, Forc, Data, v0, ode45options)

% Returns a scalar 'cost' describing model misfit to data given parameter
% updates 'x', and a struct 'costComponents' containing the cost ascribed
% to each separate data type

nt = FixedParams.nt;
nz = FixedParams.nz;
nPP_size = FixedParams.nPP_size;
nPP_nut = FixedParams.nPP_nut;
nOM_type = FixedParams.nOM_type;
nOM_nut = FixedParams.nOM_nut;
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
out.P = reshape(OUT(FixedParams.PP_index,:,:), [nPP_size nz nPP_nut nt nTraj]);
out.Z = reshape(OUT(FixedParams.ZP_index,:,:), [1 nz nt nTraj]);
out.OM = reshape(OUT(FixedParams.OM_index,:,:), [nOM_type nz nOM_nut nt nTraj]);

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
cost_POC = nan(1,nEvent);
cost_Chl = nan(1,nEvent);
for i = 1:nEvent
    iObs = Data.Event == i; % index data
    time = Data.Yearday(find(iObs,1)); % sample time    
    ti = Data.EventTraj(i,:); % index trajectories
    
    vars = unique(Data.Variable(Data.Event == i));
    
    % inorganic nitrogen
    if any(strcmp('N',vars))
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
    if any(strcmp('PON',vars))
        ind = iObs & strcmp('PON',Data.Variable);
        y_obs = Data.scaled_Value(ind); % measured values
        depths_obs = Data.Depth(ind); % measurement depths
        y_mod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_N_index ,time,ti)); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scaleFun_PON(Data.scale_mu(ind), Data.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_PON(i) = sum(sqErr(:)) / numel(y_mod);
    end    

    % POC
    if any(strcmp('POC',vars))
        ind = iObs & strcmp('POC',Data.Variable);
        y_obs = Data.scaled_Value(ind); % measured values
        depths_obs = Data.Depth(ind); % measurement depths
        y_mod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_C_index ,time,ti)); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scaleFun_POC(Data.scale_mu(ind), Data.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_POC(i) = sum(sqErr(:)) / numel(y_mod);
    end

    % Chl
    if any(strcmp('chl_a',vars))
        ind = iObs & strcmp('chl_a',Data.Variable);
        y_obs = Data.scaled_Value(ind); % measured values
        depths_obs = Data.Depth(ind); % measurement depths
        y_mod = squeeze(sum(out.P(:,:,FixedParams.PP_Chl_index,time,ti))); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scaleFun_chl_a(Data.scale_mu(ind), Data.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_Chl(i) = sum(sqErr(:)) / numel(y_mod);
    end
    
end

costComponents.N = nanmean(cost_N);
costComponents.PON = nanmean(cost_PON);
costComponents.POC = nanmean(cost_POC);
costComponents.Chl = nanmean(cost_Chl);

cost = 0;
fields = fieldnames(costComponents);
for i = 1:length(fields)
    cost = cost + costComponents.(fields{i});
end


