function [dvdt, extraOutput] = ODEs(t, v_in, parameterList)

fixedParams = parameterList.FixedParams;
params = parameterList.Params;
forc = parameterList.Forc;

%% MODEL DIMENSIONS

nz = fixedParams.nz;    % number of depth layers
nPP = fixedParams.nPP;  % number of phytoplankton size classes

%% INITIAL CONDITIONS

% Inorganic nutrient
N = v_in(fixedParams.IN_index)';

% Plankton
PP = reshape(v_in(fixedParams.PP_index), [nPP nz]);  % phytoplankton
ZP = v_in(fixedParams.ZP_index)';                    % zooplankton
B = [PP; ZP];                                        % all plankton

% Organic matter
OM = v_in(fixedParams.OM_index)';


%% FORCING DATA

% Linearly interpolate between days
T = (forc.T(:,1) + diff(forc.T,1,2) .* t)';
K = forc.K(:,1) + diff(forc.K,1,2) .* t;
I = (forc.PAR(:,1) + diff(forc.PAR,1,2) .* t)';

%% MODEL EQUATIONS

%~~~~~~~~~~~
% Physiology
%~~~~~~~~~~~

% Background mortality
B_mortality = params.m .* B;

% Temperature dependence
gammaT = exp(params.A .* (T - params.Tref));


%~~~~~~~~~~~
% Autotrophy
%~~~~~~~~~~~

% Growth depends on cellular nitrogen uptake rate and assimilation (photosynthetic) rate
if all(I == 0)
    B_uptake = zeros(nPP+1,nz);
    N_uptake_losses = zeros(1, nz);    
else    
    mu = params.pmax .* (1 - exp(-(params.aP .* I) ./ params.pmax));        % photosynthetic (metabolic) rate (1 / day)    
    mumax = params.Vmax_over_qmin ./ ((params.Qmax_over_delQ .* ... 
        params.Vmax_over_qmin) ./ mu + 1);                                  % maximum growth rate (1 / day) depends on nutrient uptake and photosynthetic rates
    V = (gammaT .* N) ./ (N ./ mumax + 1 ./ params.aN_over_qmin);           % growth rate, Michaelis-Menton form unusually written for efficiency
    B_uptake = [V; zeros(1, nz)] .* B;
    N_uptake_losses = sum(B_uptake);
end

%~~~~~~~~~~~~~
% Heterotrophy
%~~~~~~~~~~~~~

% Grazing - single predator class, with cannibalism
F = sum(B);
B2 = B.^2;
Phi = B2 ./ sum(B2);

G = params.Gmax .* gammaT .* F ./ (params.k_G + F) .* ... 
    (1 - exp(params.Lambda .* F)) .* Phi;                    % grazing rate

B_predation_losses = ZP .* G;

% Assimilation
B_predation_gains = [zeros(nPP, nz); ... 
    params.lambda_max .* sum(B_predation_losses)];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of organic matter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Only DOM is explicitly modelled as a state variable. POM is modelled to 
% sink and remineralise instantaneously.

% Messy feeding
lambda_m = 1 - params.lambda_max;
beta_x_G = params.beta .* G;
ZP_x_lambda_m = lambda_m .* ZP;
DOM_mess = ZP_x_lambda_m .* sum(beta_x_G);
POM_mess = ZP_x_lambda_m .* sum(G - beta_x_G);

% Mortality
beta_x_m_x_B = params.beta .* B_mortality;
DOM_mort = sum(beta_x_m_x_B);
POM_mort = sum(B_mortality - beta_x_m_x_B);

% Remineralisation
remin_DOM = params.rDOM .* OM;

SDOM = DOM_mort + DOM_mess - remin_DOM ; % sources/sinks of DOM
SPOM = POM_mort + POM_mess;              % sources of POM


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% POM sinking and remineralisation
remin_POM = (params.POM_to_IN * (SPOM(:) .* fixedParams.zwidth)) ./ ... 
    fixedParams.zwidth;  % to conserve mass, sinking of POM applies to quantity not concentration

% Remineralisation
SN = remin_DOM(:) + remin_POM; % with single nutrient model the only source of N is remineralised DOM and POM


%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

v_diffuse = diffusion_1D([N(:), B', OM(:)], K, ... 
    fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_diffuse = v_diffuse(:,2:nPP+2)';
OM_diffuse = v_diffuse(:,end);


%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - N_uptake_losses(:) + SN;
% Plankton
dBdt = B_diffuse + B_uptake + B_predation_gains - B_predation_losses - B_mortality;
% Organic matter
dOMdt = OM_diffuse + SDOM(:);

% Separate phytoplankton and zooplankton
dPPdt = dBdt(fixedParams.phytoplankton,:);
dZPdt = dBdt(fixedParams.zooplankton,:);

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt];





%% AUXILIARY OUTPUTS
if fixedParams.returnExtras
    extraOutput.POM = SPOM;
    extraOutput.remin_POM = remin_POM';
end

end


%% Functions

function v = diffusion_1D(u,K,w,delz)
% Returns rate of change of u due to diffusion.
% Inputs: array of concentrations u,
%         array of diffusivities K,
%         depth layer widths w,
%         distance between depth layer centers delz.
%         The 1st dimension of each array is space (modelled depths of
%         water column).
% Function assumes zero flux boundary conditions.
% Diffusivities K are at midpoints of spatial grid used for u, so 
% size(K,1)=size(u,1)-1.
usize = size(u);
dmdt = zeros(usize);
dmdt(1:usize(1)-1,:) = diff([zeros(1,usize(2)); (K ./ delz) .* diff(u)]); % rate of change of mass
dmdt(usize(1),:) = -sum(dmdt);  % mass is conserved
v = dmdt ./ w; % rate of change of concentration
end

