function v0 = initialiseVariables(FixedParams, Params, Forc, varargin)
% Set initial values, v0, for all state variables.
% Uses SINMOD estimates of NO3 and planktonic N concentrations and some
% assumptions to derive initial values. Assumptions involve using some
% model parameter values.

extractVarargin(varargin)

v0 = nan(FixedParams.nEquations, Forc.nTraj);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Inorganic nutrients - dimension = [depth, trajectory]

NO3ic = squeeze(Forc.NO3ic(:,1,:));
v0(FixedParams.IN_index,:) = NO3ic;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Plankton - dimension = [size, depth, nutrient, trajectory]

% Split the SINMOD (PS = small, PL = large) plankton N estimates evenly
% across modelled size classes
smallCell_N = squeeze(Forc.PSic(:,1,:)); % dimension: [depth, trajectory]
largeCell_N = squeeze(Forc.PLic(:,1,:));
PS = ~FixedParams.diatoms; % index small & large cells
PL = FixedParams.diatoms;
nSmall = sum(PS);
nLarge = sum(PL);
% Split initial plankton N evenly across modelled size classes -- grouped
% by small or large
smallCell_N = (1/nSmall) .* repmat(reshape(smallCell_N, ...
    [1 size(smallCell_N)]), [nSmall 1 1]);
largeCell_N = (1/nLarge) .* repmat(reshape(largeCell_N, ...
    [1 size(largeCell_N)]), [nLarge 1 1]); % dimension: [size, depth, trajectory]

% Indexes matching SINMOD outputs to modelled variables
sp = 1:sum(PS(FixedParams.phytoplankton));
lp = 1:sum(PL(FixedParams.phytoplankton));
sz = sum(PS(FixedParams.phytoplankton))+1:size(smallCell_N, 1);
lz = sum(PL(FixedParams.phytoplankton))+1:size(largeCell_N, 1);

% Store all plankton initial values in struct P.
P.N = cat(1, smallCell_N(sp,:,:), largeCell_N(lp,:,:), ...
    smallCell_N(sz,:,:), largeCell_N(lz,:,:)); % merge arrays of small and large plankton

% Initialise carbon by assuming N quotas are 100*Nquota% full (0 <= Nquota <= 1)
if exist('Nquota', 'var') && (eval('Nquota') < 0 || eval('Nquota') > 1)
    Nquota = 0.75;
    warning(['Optional argument "Nquota" out of bounds: constraint is 0<=Nquota<=1. Value reset to default Nquota=' num2str(Nquota)])
end
if ~exist('Nquota', 'var')
    Nquota = 0.75; % N quota at 3/4 reasonable for initial date of Jan 1st
end

Q_N = Params.Qmin_QC + Nquota .* Params.delQ_QC;
P.C = P.N ./ Q_N;

% Initialise chlorophyll as chlFrac the maximum ratio with N
if exist('chlFrac', 'var') && (eval('chlFrac') < 0 || eval('chlFrac') > 1)
    chlFrac = 0.75;
    warning(['Optional argument "chlFrac" out of bounds: constraint is 0<=chlFrac<=1. Value reset to default chlFrac=' num2str(chlFrac)])
end
if ~exist('chlFrac', 'var')
    % Chl/N ratio at 3/4 the maximum reasonable for initial date of Jan 1st    
    chlFrac = 0.75;
end

chlN_max = Params.theta; % maximum Chl:N ratio (mg Chl / mmol N)
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

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Organic matter - dimension = [type, depth, nutrient, trajectory]

OM = nan(FixedParams.nOM_type, FixedParams.nz, FixedParams.nOM_nut, Forc.nTraj);

% No SINMOD output for organic matter => initialise by assumed ratios with
% plankton. Use OM_frac>0 to specify quantity of OM relative to plankton,
% and use 0<DOM_frac<1 to specify proportion of DOM relative to POM.
if exist('OM_frac', 'var') && eval('OM_frac') < 0
    OM_frac = 0.05;
    warning(['Optional argument "OM_frac" out of bounds: constraint is 0<OM_frac. Value reset to default OM_frac=' num2str(OM_frac)])
end
if ~exist('OM_frac', 'var')
    OM_frac = 0.05;
end
if exist('DOM_frac', 'var') && eval('DOM_frac') < 0
    DOM_frac = 0.5;
    warning(['Optional argument "DOM_frac" out of bounds: constraint is 0<DOM_frac<1. Value reset to default DOM_frac=' num2str(DOM_frac)])
end
if ~exist('DOM_frac', 'var')
    DOM_frac = 0.5;
end

tot_PN = squeeze(sum(P.N));
tot_PC = squeeze(sum(P.C));

OC = OM_frac .* tot_PC;
DOC = DOM_frac * OC;
POC = OC - DOC;
ON = OM_frac .* tot_PN;
DON = DOM_frac * ON;
PON = ON - DON;
OM(FixedParams.DOM_index,:,FixedParams.OM_C_index,:) = DOC;
OM(FixedParams.DOM_index,:,FixedParams.OM_N_index,:) = DON;
OM(FixedParams.POM_index,:,FixedParams.OM_C_index,:) = POC;
OM(FixedParams.POM_index,:,FixedParams.OM_N_index,:) = PON;

v0(FixedParams.OM_index,:) = reshape(OM, [FixedParams.nOM * FixedParams.nz Forc.nTraj]);
