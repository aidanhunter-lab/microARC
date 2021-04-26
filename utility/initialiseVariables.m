function v0 = initialiseVariables(FixedParams, Params, Forc)

v0 = nan(FixedParams.nEquations, Forc.nTraj); % state variable initial values
% Use SINMOD forcing to set initial values

% Inorganic nutrients - dimension = [depth, trajectory]
NO3ic = squeeze(Forc.NO3ic(:,1,:));
v0(FixedParams.IN_index,:) = NO3ic;

% Plankton - dimension = [size, depth, nutrient, trajectory]
% Split the SINMOD (PS & PL) plankton N estimates evenly across length classes
smallCell_N = squeeze(Forc.PSic(:,1,:));
largeCell_N = squeeze(Forc.PLic(:,1,:));
nSmall = sum(~FixedParams.diatoms);
nLarge = sum(FixedParams.diatoms);
smallCell_N = (1/nSmall) .* repmat(reshape(smallCell_N, ...
    [1 size(smallCell_N)]), [nSmall 1 1]);
largeCell_N = (1/nLarge) .* repmat(reshape(largeCell_N, ...
    [1 size(largeCell_N)]), [nLarge 1 1]);

PS = FixedParams.phytoplankton & ~FixedParams.diatoms;
PL = FixedParams.phytoplankton & FixedParams.diatoms;

P.N = cat(1, smallCell_N(1:sum(PS),:,:), largeCell_N(1:sum(PL),:,:), ...
    smallCell_N(sum(PS)+1:end,:,:), largeCell_N(sum(PL)+1:end,:,:)); % merge arrays of small and large plankton


% Initialise carbon by assuming N quotas are 0<Crat<1 times the maximum
Crat = 0.75;
Q_N = Params.Qmin_QC + Crat .* Params.delQ_QC;
P.C = P.N ./ Q_N;

% Set initial chlorophyll as chlFrac the maximum ratio with N
chlN_max = Params.theta; % maximum Chl:N ratio (mg Chl / mmol N)
chlFrac = 0.75;
% P.Chl = chlFrac * chlN_max .* P.N(FixedParams.phytoplankton,:,:);
P.Chl = chlFrac * chlN_max .* P.N;
P.Chl(FixedParams.zooplankton,:,:) = nan;

% autotrophs
PP = structfun(@(z) z(FixedParams.phytoplankton,:,:), P, 'UniformOutput', false);
vPP = nan(FixedParams.nPP_size, FixedParams.nz, FixedParams.nPP_nut, Forc.nTraj);
for i = 1:FixedParams.nPP_nut
    vPP(:,:,i,:) = reshape(PP.(FixedParams.PP_nut{i}), ...
        [FixedParams.nPP_size FixedParams.nz 1 Forc.nTraj]);
end

v0(FixedParams.PP_index,:) = reshape(vPP, [FixedParams.nPP * FixedParams.nz Forc.nTraj]);

% heterotrophs
ZP = structfun(@(z) z(FixedParams.zooplankton,:,:), P, 'UniformOutput', false);
vZP = nan(FixedParams.nZP_size, FixedParams.nz, FixedParams.nZP_nut, Forc.nTraj);
for i = 1:FixedParams.nZP_nut
    vZP(:,:,i,:) = reshape(ZP.(FixedParams.ZP_nut{i}), ...
        [FixedParams.nZP_size FixedParams.nz 1 Forc.nTraj]);
end

v0(FixedParams.ZP_index,:) = reshape(vZP, [FixedParams.nZP * FixedParams.nz Forc.nTraj]);


% Organic matter - dimension = [type, depth, nutrient, trajectory]
OM = nan(FixedParams.nOM_type, FixedParams.nz, FixedParams.nOM_nut, Forc.nTraj);
tot_PN = squeeze(sum(P.N));
tot_PC = squeeze(sum(P.C));

OM_frac = 0.05;
DOM_frac = 0.5;
DOC = DOM_frac * OM_frac .* tot_PC;
POC = (1-DOM_frac) * OM_frac .* tot_PC;
DON = DOM_frac * OM_frac .* tot_PN;
PON = (1-DOM_frac) * OM_frac .* tot_PN;
OM(FixedParams.DOM_index,:,FixedParams.OM_C_index,:) = DOC;
OM(FixedParams.DOM_index,:,FixedParams.OM_N_index,:) = DON;
OM(FixedParams.POM_index,:,FixedParams.OM_C_index,:) = POC;
OM(FixedParams.POM_index,:,FixedParams.OM_N_index,:) = PON;

v0(FixedParams.OM_index,:) = reshape(OM, [FixedParams.nOM * FixedParams.nz Forc.nTraj]);


