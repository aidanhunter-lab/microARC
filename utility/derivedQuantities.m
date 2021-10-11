function auxVars = derivedQuantities(Forc, FixedParams, Params, out, auxVars)

dt = unique(diff(Forc.t(:,1)));
nTraj = Forc.nTraj;

extractStruct(FixedParams)
nsize = nPP_size + nZP_size;

% Primary production
auxVars.PP_uptake = out.P(:,:,PP_C_index,:,:) .* ...
    auxVars.V(phytoplankton,:,:,:,:); % mmol x / m3 / day
auxVars.PP_uptake_depthInt = sum(reshape(zwidth, [1 nz]) .* auxVars.PP_uptake , 2); % depth-integrated production, mmol x / m2 / day
auxVars.PP_uptake_tot = sum(dt .* auxVars.PP_uptake_depthInt, ...
    find(size(auxVars.PP_uptake_depthInt) == nt)); % total production of each size class along each trajectory, mmol C / m2 / year (simulation time)

% Feeding fluxes (grazing loss) and secondary production (grazing gain)
Z_C = out.Z(:,:,ZP_C_index,:,:); % zooplankton carbon
auxVars.predation_losses = reshape(auxVars.Q(:,:,~PP_Chl_index,:,:), ...
    [1, nPP_size + nZP_size, nz, nZP_nut, nt, nTraj]) .* ...
    reshape(Z_C, [nZP_size, 1, nz, 1, nt, nTraj]) .* ...
    reshape(auxVars.G, [nZP_size, nPP_size + nZP_size, nz, 1, nt, nTraj]); % mmol x / m3 / day
auxVars.predation_losses_depthInt = sum(reshape(zwidth, [1, 1, nz]) .* ...
    auxVars.predation_losses, 3, 'omitnan'); % mmol x / m2 / day
auxVars.predation_losses_tot = sum(dt .* auxVars.predation_losses_depthInt, ...
    find(size(auxVars.predation_losses_depthInt) == nt), 'omitnan'); % mmol x / m2 / year

predation_gains_all = auxVars.lambda(:,:,:,~PP_Chl_index,:,:) .* auxVars.predation_losses;
mess = auxVars.predation_losses - predation_gains_all;

auxVars.predation_gains = sum(predation_gains_all, 2, 'omitnan');
auxVars.predation_gains_depthInt = sum(reshape(zwidth, [1, 1, nz]) .* ...
    auxVars.predation_gains, 3, 'omitnan');
auxVars.predation_gains_tot = sum(dt .* auxVars.predation_gains_depthInt, ...
    find(size(auxVars.predation_gains_depthInt) == nt), 'omitnan');

% OM production
% Messy feeding
beta_mess = reshape(Params.beta, [1, nsize]) .* mess;
auxVars.OM_mess_DOM = sum(beta_mess, 2);
auxVars.OM_mess_POM = sum(mess - beta_mess, 2);
auxVars.OM_mess_DOM_depthInt = sum(reshape(zwidth, [1, 1, nz]) .* ...
    auxVars.OM_mess_DOM, 3);
auxVars.OM_mess_DOM_tot = sum(dt .* auxVars.OM_mess_DOM_depthInt, ...
    find(size(auxVars.OM_mess_DOM_depthInt) == nt));
auxVars.OM_mess_POM_depthInt = sum(reshape(zwidth, [1, 1, nz]) .* ...
    auxVars.OM_mess_POM, 3);
auxVars.OM_mess_POM_tot = sum(dt .* auxVars.OM_mess_POM_depthInt, ...
    find(size(auxVars.OM_mess_POM_depthInt) == nt));
% Mortality
B = [out.P(:,:,~PP_Chl_index,:,:); out.Z];
mortality = (Params.m + Params.m2 .* B) .* B;
auxVars.OM_mort_DOM = Params.beta .* mortality;
auxVars.OM_mort_POM = mortality - auxVars.OM_mort_DOM;
auxVars.OM_mort_DOM_depthInt = sum(reshape(zwidth, [1, nz]) .* ...
    auxVars.OM_mort_DOM, 2);
auxVars.OM_mort_DOM_tot = sum(dt .* auxVars.OM_mort_DOM_depthInt, ...
    find(size(auxVars.OM_mort_DOM_depthInt) == nt));
auxVars.OM_mort_POM_depthInt = sum(reshape(zwidth, [1, nz]) .* ...
    auxVars.OM_mort_POM, 2);
auxVars.OM_mort_POM_tot = sum(dt .* auxVars.OM_mort_POM_depthInt, ...
    find(size(auxVars.OM_mort_POM_depthInt) == nt));


