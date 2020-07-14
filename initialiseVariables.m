function v0 = initialiseVariables(FixedParams, Params, Forc)

v0 = nan(FixedParams.nEquations, Forc.nTraj); % state variable initial values
% Use SINMOD forcing to set initial values

% Inorganic nutrients - dimension = [depth, trajectory]
NO3ic = squeeze(Forc.NO3ic(:,1,:));
v0(FixedParams.IN_index,:) = NO3ic;

% Phytoplankton - dimension = [size, depth, nutrient, trajectory]
% Split the SINMOD (PS & PL) plankton N estimates evenly across length classes
smallCell_N = squeeze(Forc.PSic(:,1,:));
largeCell_N = squeeze(Forc.PLic(:,1,:));
nSmall = sum(~FixedParams.diatoms);
nLarge = sum(FixedParams.diatoms);
smallCell_N = (1/nSmall) .* repmat(reshape(smallCell_N, ...
    [1 size(smallCell_N)]), [nSmall 1 1]);
largeCell_N = (1/nLarge) .* repmat(reshape(largeCell_N, ...
    [1 size(largeCell_N)]), [nLarge 1 1]);
P.N = cat(1,smallCell_N,largeCell_N);                        % merge arrays of small and large plankton

% Initialise carbon by assuming N quotas are 0<Crat<1 the maximum
Crat = 0.75;
Q_N = Params.Qmin_QC + Crat .* Params.delQ_QC;
P.C = P.N ./ Q_N;

% % Initialise carbon by assuming a C:N ratio of 163:22 (Martiny et al. Sci. Data 1, 140048 (2014))
% CNrat = 163/22;
% P.C = CNrat * P.N;

% Set initial chlorophyll as chlFrac the maximum ratio with N
chlN_max = Params.theta; % maximum Chl:N ratio (mg Chl / mmol N)
chlFrac = 0.75;
P.Chl = chlFrac * chlN_max .* P.N;

% combine C, N & Chl
P0 = nan(FixedParams.nPP_size, FixedParams.nz, FixedParams.nPP_nut, Forc.nTraj);
for i = 1:FixedParams.nPP_nut
    P0(:,:,i,:) = reshape(P.(FixedParams.PP_nut{i}), ...
        [FixedParams.nPP_size FixedParams.nz 1 Forc.nTraj]);
%     P0(:,:,i,:) = reshape(eval(['P_' FixedParams.PP_nut{i}]), ... 
%         [FixedParams.nPP_size FixedParams.nz 1 Forc.nTraj]);
end

v0(FixedParams.PP_index,:) = reshape(P0, [FixedParams.nPP * FixedParams.nz Forc.nTraj]);


% Zooplankton - dimension = [depth, trajectory]
% Zooplankton are not assigned any particular size class, so just
% initialise their total carbon content as some fraction of the total
% phytoplankton carbon content
Z_P_frac = 0.25;
tot_PC = squeeze(sum(P.C));
v0(FixedParams.ZP_index,:) = Z_P_frac * tot_PC;

% Organic matter - dimension = [type, depth, nutrient, trajectory]
OM = nan(FixedParams.nOM_type, FixedParams.nz, FixedParams.nOM_nut, Forc.nTraj);
tot_PN = squeeze(sum(P.N));

OM_frac = 0.05;
DOM_frac = 0.5;
DOC = DOM_frac * OM_frac .* (v0(FixedParams.ZP_index,:) + tot_PC);
POC = (1-DOM_frac) * OM_frac .* (v0(FixedParams.ZP_index,:) + tot_PC);
DON = DOM_frac * OM_frac .* (Z_P_frac * tot_PN + tot_PN);
PON = (1-DOM_frac) * OM_frac .* (Z_P_frac * tot_PN + tot_PN);
OM(FixedParams.DOM_index,:,FixedParams.OM_C_index,:) = DOC;
OM(FixedParams.DOM_index,:,FixedParams.OM_N_index,:) = DON;
OM(FixedParams.POM_index,:,FixedParams.OM_C_index,:) = POC;
OM(FixedParams.POM_index,:,FixedParams.OM_N_index,:) = PON;

v0(FixedParams.OM_index,:) = reshape(OM, [FixedParams.nOM * FixedParams.nz Forc.nTraj]);


