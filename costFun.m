function [cost, costComponents] = costFun(x, FixedParams, Params, Forc, Data, v0, odeIntegrator, odeOptions)

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
nEvent = Data.scalar.nEvents;
depths_mod = abs(FixedParams.z);

% Set parameter values
Params = updateParameters(Params, FixedParams, x);

% Integrate
[OUT, AUXVARS, AUXVARS_2d, namesExtra, nExtra] = ... 
    integrateTrajectories(FixedParams, Params, Forc, v0, odeIntegrator, odeOptions);

OUT = exp(OUT); % exponentiate to natural scale

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
    iObs = Data.scalar.Event == i; % index data
    time = Data.scalar.Yearday(find(iObs,1)); % sample time    
    ti = Data.scalar.EventTraj(i,:); % index trajectories
    
    vars = unique(Data.scalar.Variable(Data.scalar.Event == i));
    
    % inorganic nitrogen
    if any(strcmp('N',vars))
        ind = iObs & strcmp('N',Data.scalar.Variable);
        y_obs = Data.scalar.scaled_Value(ind); % measured values
        depths_obs = Data.scalar.Depth(ind); % measurement depths
        y_mod = squeeze(out.N(:,:,time,ti)); % modelled values [depth,trajectory]
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths        
        y_mod = Data.scalar.scaleFun_N(Data.scalar.scale_mu(ind), ...
            Data.scalar.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_N(i) = sum(sqErr(:)) / numel(y_mod);
    end
    
    % PON
    if any(strcmp('PON',vars))
        ind = iObs & strcmp('PON',Data.scalar.Variable);
        y_obs = Data.scalar.scaled_Value(ind); % measured values
        depths_obs = Data.scalar.Depth(ind); % measurement depths
        y_mod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_N_index ,time,ti)); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scalar.scaleFun_PON(Data.scalar.scale_mu(ind), Data.scalar.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_PON(i) = sum(sqErr(:)) / numel(y_mod);
    end    

    % POC
    if any(strcmp('POC',vars))
        ind = iObs & strcmp('POC',Data.scalar.Variable);
        y_obs = Data.scalar.scaled_Value(ind); % measured values
        depths_obs = Data.scalar.Depth(ind); % measurement depths
        y_mod = squeeze(out.OM(FixedParams.POM_index,:,FixedParams.OM_C_index ,time,ti)); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scalar.scaleFun_POC(Data.scalar.scale_mu(ind), Data.scalar.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_POC(i) = sum(sqErr(:)) / numel(y_mod);
    end

    % Chl
    if any(strcmp('chl_a',vars))
        ind = iObs & strcmp('chl_a',Data.scalar.Variable);
        y_obs = Data.scalar.scaled_Value(ind); % measured values
        depths_obs = Data.scalar.Depth(ind); % measurement depths
        y_mod = squeeze(sum(out.P(:,:,FixedParams.PP_Chl_index,time,ti))); % modelled values
        y_mod = interp1(depths_mod, y_mod, depths_obs); % interpolate modelled output to match observation depths
        y_mod = Data.scalar.scaleFun_chl_a(Data.scalar.scale_mu(ind), Data.scalar.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_Chl(i) = sum(sqErr(:)) / numel(y_mod);
    end
    
end

costComponents.N = nanmean(cost_N);
costComponents.PON = nanmean(cost_PON);
costComponents.POC = nanmean(cost_POC);
costComponents.Chl = nanmean(cost_Chl);

% Size spectra cost component
cost_N_at_size = nan(1, nPP_size);
% esd = unique(Data.size.size);
for i = 1:nPP_size
    
    iObs = Data.size.sizeClass == i; % index data
    time = Data.size.Yearday(find(iObs,1)); % sample time

    vars = unique(Data.size.Variable);
    
    % planktonic nitrogen at size
    if any(strcmp('N_at_size',vars))
        ind = iObs & strcmp('N_at_size',Data.size.Variable);
        y_obs = Data.size.scaled_Value(ind); % measured values
        depths_obs = [Data.size.DepthMin(ind) Data.size.DepthMax(ind)]; % measurement depth range        
        depth_ind = -FixedParams.zw(2:end) > depths_obs(1,1) & ...
            -FixedParams.zw(1:end-1) < depths_obs(1,2); % depth layers corresponding to samples
        
        y_mod = squeeze(out.P(i,depth_ind,FixedParams.PP_N_index,time,:)); % modelled values
        y_mod = mean(y_mod); % average over depth layers        
        y_mod = Data.size.scaleFun(Data.size.scale_mu(ind), ...
            Data.size.scale_sig(ind), y_mod); % scale model output using same functions that scaled the data
        sqErr = (y_obs - y_mod) .^2; % normal error
        cost_N_at_size(i) = sum(sqErr(:)) / numel(sqErr);
    end
end

costComponents.N_at_size = nanmean(cost_N_at_size);

cost = sum(struct2array(costComponents));

% cost = 0;
% fields = fieldnames(costComponents);
% for i = 1:length(fields)
%     cost = cost + costComponents.(fields{i});
% end


