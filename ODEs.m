function [dvdt, extraOutput, extraOutput_2d] = ODEs(t, v_in, parameterList, forc, timeStep, returnExtra)

fixedParams = parameterList.FixedParams;
params = parameterList.Params;

%% MODEL DIMENSIONS

nz = fixedParams.nz;    % number of depth layers
nPP_size = fixedParams.nPP_size;  % number of phytoplankton size classes
nPP_nut = fixedParams.nPP_nut;  % number of phytoplankton nutrient classes
nOM_type = fixedParams.nOM_type;  % number of organic matter types
nOM_nut = fixedParams.nOM_nut;  % number of organic nutrient classes
phyto = fixedParams.phytoplankton;
zoo = fixedParams.zooplankton;

%% INITIAL CONDITIONS

% v_in_exp = exp(v_in); % exponentiate to natural scale

% Inorganic nitrogen
% N = v_in_exp(fixedParams.IN_index)';
N = v_in(fixedParams.IN_index)';

% Plankton
% PP = reshape(v_in_exp(fixedParams.PP_index), [nPP_size nz nPP_nut]); % phytoplankton (all nutrients)
PP = reshape(v_in(fixedParams.PP_index), [nPP_size nz nPP_nut]); % phytoplankton (all nutrients)
P_C = PP(:,:,fixedParams.PP_C_index);
% Z_C = v_in_exp(fixedParams.ZP_index)'; % zooplankton (carbon)
Z_C = v_in(fixedParams.ZP_index)'; % zooplankton (carbon)
B_C = [P_C; Z_C]; % all planktonic carbon

% Organic matter
% OM =reshape(v_in_exp(fixedParams.OM_index), [nOM_type nz nOM_nut]);
OM =reshape(v_in(fixedParams.OM_index), [nOM_type nz nOM_nut]);
OM_C = OM(:,:,fixedParams.OM_C_index); % DOC and POC


%% FORCING DATA

T = forc.T(:,timeStep-1:timeStep);
K = forc.K(:,timeStep-1:timeStep);
Isurf = forc.PARsurf(:,timeStep-1:timeStep);

% Linearly interpolate between days
T = (T(:,1) + diff(T,1,2) .* t)';
K = K(:,1) + diff(K,1,2) .* t;
Isurf = (Isurf(:,1) + diff(Isurf,1,2) .* t)';

% Calculate light levels at depth -- within each depth layer light
% attenuates over half the layer width plus the combined widths of all
% shallower layers.
att = (fixedParams.attSW + fixedParams.attP .* sum(PP(:,:,fixedParams.PP_Chl_index))) .* fixedParams.zwidth';
att = 0.5 * att + [0 cumsum(att(1:nz-1))];
I = Isurf * exp(-att);


%% MODEL EQUATIONS

%~~~~~~~~~~~
% Physiology
%~~~~~~~~~~~

% Phytoplankton N and Chl quotas relative to C
Q_N = PP(:,:,fixedParams.PP_N_index) ./ P_C;
Q_Chl = PP(:,:,fixedParams.PP_Chl_index) ./ P_C;
% Q_Chl_N = PP(:,:,fixedParams.PP_Chl_index) ./ PP(:,:,fixedParams.PP_N_index); % Chl:N ratio

% Nutrient limitation
gammaN = max(0, min(1, (Q_N - params.Qmin_QC) ./ params.delQ_QC));

% Uptake regulation
Qstat = 1 - gammaN .^ params.h;

% Temperature dependence
gammaT = exp(params.A .* (T - params.Tref));

% Background mortality
B_C_mortality = params.m .* B_C; % linear mortality
% B_C_mortality = params.m .* B_C .^ 2; % non-linear mortality


%~~~~~~~~~~~
% Autotrophy
%~~~~~~~~~~~

zeros_size_nz = zeros(nPP_size, nz);

% Nutrient uptake
V_N = MichaelisMenton(params.Vmax_QC, params.kN, N) .* gammaT .* Qstat;

% V_N = params.Vmax_QC ./ (1 + params.Vmax_QC ./ (params.aN_QC .* N)) .* gammaT .* Qstat; % nitrogen uptake rate (mmol N / mmol C / day)
N_uptake = V_N .* P_C; % mmol N / m^3 / day
N_uptake_losses = sum(N_uptake);

% Photosynthesis
if all(I == 0)
    V_Chl = zeros_size_nz;
    V_C = zeros(nPP_size+1, nz);
else
    psat = params.pmax .* gammaT .* gammaN; % light saturated photosynthetic rate
    aP_Q_I = (params.aP .* I) .* Q_Chl;
    pc = psat .* (1 - exp(-aP_Q_I ./ psat )); % photosynthetic (carbon production) rate (1 / day)
    rho = params.theta .* pc ./ aP_Q_I;  % proportion of new nitrogen prodcution allocated to chlorophyll (mg Chl / mmol N)
    V_Chl = rho .* V_N; % chlorophyll production rate (mg Chl / mmol C / day)
    V_C = [max(0, pc - params.xi .* V_N); zeros(1, nz)];
end


%~~~~~~~~~~~~~
% Heterotrophy
%~~~~~~~~~~~~~

% Grazing - single predator class, with cannibalism
F = sum(B_C); % total prey carbon
BC2 = B_C .^ 2;
Phi = BC2 ./ sum(BC2); % prey preference
G = (MichaelisMenton(params.Gmax, params.k_G, F) .* ... 
    gammaT .* (1-exp(params.Lambda .* F))) .* Phi;  % grazing rate (1 / day)

% G = (params.Gmax .* gammaT .* F ./ (params.k_G + F) .* ... 
%     (1-exp(params.Lambda .* F))) .* Phi; % grazing rate (1 / day)

predation_losses_C = G .* Z_C; % mmol C / m^3 / day
predation_gains_C = [zeros_size_nz; ... 
    params.lambda_max .* sum(predation_losses_C)];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of organic matter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% OM sources are phyto and zooplankton mortality. Organic nitrogen is
% modelled but there is no zooplankton N class, so assume that zooplankton
% N quota is the same as their prey (some data on ths is would be useful...)
phyto_C_losses = predation_losses_C(phyto,:); % Assume that zooplankton N quota is the same as their prey (some data on ths is would be useful...)
Q_N_zoo = sum(Q_N .* (phyto_C_losses ./ sum(phyto_C_losses)));

% Messy feeding
lambda_predLoss = (1 - params.lambda_max) .* predation_losses_C;
beta_lambda_predLoss = params.beta .* lambda_predLoss;
OM_mess_C = [sum(beta_lambda_predLoss); ... 
    sum(lambda_predLoss-beta_lambda_predLoss)]; % OM_mess=[DOM_mess; POM_mess] (mmol C / m^3 / day)
OM_mess_N = Q_N_zoo .* OM_mess_C;

% Mortality
beta_m_B = params.beta .* B_C_mortality;
OM_mort_C = [sum(beta_m_B); sum(B_C_mortality - beta_m_B)]; % OM_mort=[DOM_mort; POM_mort] (mmol C / m^3 / day)
B_N_mortality = [Q_N; Q_N_zoo] .* B_C_mortality;
beta_m_B = params.beta .* B_N_mortality;
OM_mort_N = [sum(beta_m_B); sum(B_N_mortality - beta_m_B)]; % OM_mort=[DOM_mort; POM_mort] (mmol N / m^3 / day)

% Remineralisation
OM_remin = params.rOM .* OM;
OM_remin_N = OM_remin(:,:,fixedParams.OM_N_index);

SOM_C = OM_mort_C + OM_mess_C - OM_remin(:,:,fixedParams.OM_C_index); % (mmol C / m^3 / day)
SOM_N = OM_mort_N + OM_mess_N - OM_remin_N; % (mmol N / m^3 / day)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Remineralisation
SN = sum(OM_remin_N); % (mmol N / m^3 / day)

%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

B_C_t = B_C';
OM_C_t = OM_C';

v_diffuse = diffusion_1D([N(:), B_C_t, OM_C_t], K, fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_C_diffuse = v_diffuse(:,2:nPP_size+2)';
OM_C_diffuse = v_diffuse(:,end-nOM_nut+1:end)';

%~~~~~~~~
% Sinking
%~~~~~~~~

v_sink = sinking([B_C_t(:,phyto), OM_C_t], [params.wp, params.wk], fixedParams.zwidth);

B_C_sink = [v_sink(:,1:nPP_size), zeros(nz, 1)]';
OM_C_sink = v_sink(:,nPP_size+1:end)';

%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - N_uptake_losses(:) + SN(:);
% Plankton
fluxC_ = B_C_sink + B_C_diffuse - predation_losses_C - B_C_mortality; % terms applicable to both phyto and zoo - all plankton carbon flux terms (mmol C / m^3 / day) excluding uptake and predation gains
fluxC = fluxC_ + predation_gains_C; % all plankton carbon flux terms excluding uptake
fluxC_p = fluxC_(phyto,:);
fluxN = N_uptake + Q_N .* fluxC_p; % mmol N / m^3 / day
fluxChl = V_Chl .* P_C + Q_Chl .* fluxC_p; % mg Chl / m^3 / day
fluxC = V_C .* B_C + fluxC;
dPPdt = cat(3, fluxC(phyto,:), fluxN, fluxChl);
dZPdt = fluxC(zoo,:);
% Organic matter
fluxC_ = OM_C_diffuse + OM_C_sink;
fluxN = OM(:,:,fixedParams.OM_N_index) ./ OM_C .* fluxC_ + SOM_N;
fluxC = fluxC_ + SOM_C;
dOMdt = cat(3, fluxC, fluxN);

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt(:)];



%% AUXILIARY OUTPUTS
if returnExtra
    % 1D
    extraOutput.PAR = I;
    % 2D (vectorised over cell size)
    extraOutput_2d = struct();
    extraOutput_2d.cellDensity = P_C ./ params.Q_C; % phytoplankton cell density [cells / m^3]
    extraOutput_2d.biovolume = (1e-18 * fixedParams.PPsize) .* ... 
        extraOutput_2d.cellDensity; % [m^3 / m^3] volume of all cells per m^3 of water
%     if all(I == 0)
%         extraOutput_2d.Qeq = repmat(params.Qmax, [1 nz]); % equilibrium nitrogen Qeq=Qmax when I=0
%     else
%         Q_ = (QmaxVmax_ .* aN) ./ (aN + params.Vmax_over_Qmin);
%         extraOutput_2d.Qeq = (Q_ + mu) ./ (mu ./ params.Qmin + Q_ ./ params.Qmax); % equilibrium N quota
%     end
%     extraOutput_2d.PP_C = params.Q_C  .* (PP(:,:,fixedParams.N_index) ./ extraOutput_2d.Qeq); % phytoplankton carbon biomass (mmol C / m^3)
end


end


%% Functions

function v = diffusion_1D(u,K,w,delz)
% Rate of change of u due to diffusion, assuming zero flux boundary conditions.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         K = diffusivities, size(K)=[nz-1 1]
%         w = depth layer widths, size(w)=[nz 1]
%         delz = distance between depth layer centers, size(delz)=[nz-1 1]
z = zeros(1,size(u,2));
v = diff([z; (K ./ delz) .* diff(u); z]) ./ w;
end

function v = sinking(u,s,w)
% Rate of change of u due to sinking.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         s = sinking speed, size(s)=[1 nvar]
%         w = depth layer widths, size(w)=[nz 1]
v = ((-s) ./ w) .* diff([zeros(1,size(u,2)); u]);
end

function v = MichaelisMenton(m,k,u)
% Uptake rate of u, given maximum m and half saturation k
u(u<0) = 0; % include for robustness... there shouldn't be any negatives
v = m .* u ./ (u + k);
end


