function v0 = initialiseVariables(FixedParams, Forc)

v0 = nan(FixedParams.nEquations, Forc.nTraj); % state variable initial values
% Use SINMOD forcing to set initial values

% Inorganic nutrients - dimension = [depth, trajectory]
NO3ic = squeeze(Forc.NO3ic(:,1,:));
v0(FixedParams.IN_index,:) = NO3ic;

% Phytoplankton - dimension = [size, depth, trajectory]
% Split the SINMOD (PS & PL) plankton estimates evenly across length classes
smallCell_N = squeeze(Forc.PSic(:,1,:));
largeCell_N = squeeze(Forc.PLic(:,1,:));
nSmall = sum(~FixedParams.diatoms);
nLarge = sum(FixedParams.diatoms);
smallCell_N = (1/nSmall) .* repmat(reshape(smallCell_N, ...
    [1 size(smallCell_N)]), [nSmall 1 1]);
largeCell_N = (1/nLarge) .* repmat(reshape(largeCell_N, ...
    [1 size(largeCell_N)]), [nLarge 1 1]);
P_N = cat(1,smallCell_N,largeCell_N);                        % merge arrays of small and large plankton
v0(FixedParams.PP_index,:) = reshape(P_N, [FixedParams.nPP * FixedParams.nz Forc.nTraj]);

% Zooplankton - dimension = [depth, trajectory]
% Zooplankton are not assigned any particular size class, so just
% initialise their total nitrogen content as some fraction of the total
% phytoplankton nitrogen content
Z_P_frac = 0.25;
tot_P = squeeze(sum(P_N));
v0(FixedParams.ZP_index,:) = Z_P_frac * tot_P;

% Organic matter - dimension = [type, depth, trajectory]
OM_frac = 0.05;
DOM_frac = 0.5;
OM = OM_frac * (v0(FixedParams.ZP_index,:) + tot_P);
DOM = DOM_frac * OM;
POM = OM - DOM;
OM = nan(FixedParams.nOM, FixedParams.nz, Forc.nTraj);
OM(FixedParams.DOM_index,:,:) = DOM;
OM(FixedParams.POM_index,:,:) = POM;
v0(FixedParams.OM_index,:) = reshape(OM, [FixedParams.nOM * FixedParams.nz Forc.nTraj]);

