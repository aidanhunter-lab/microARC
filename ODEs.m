function [dvdt, extraOutput, extraOutput_2d] = ODEs(t, v_in, parameterList, forc, timeStep, returnExtra)

fixedParams = parameterList.FixedParams;
params = parameterList.Params;

%% MODEL DIMENSIONS

nz = fixedParams.nz;    % number of depth layers
nPP = fixedParams.nPP;  % number of phytoplankton size classes
nOM = fixedParams.nOM;  % number of organic matter types

%% INITIAL CONDITIONS

% Inorganic nitrogen
N = v_in(fixedParams.IN_index)';

% Plankton
PP = reshape(v_in(fixedParams.PP_index), [nPP nz]);  % phytoplankton
ZP = v_in(fixedParams.ZP_index)';                    % zooplankton
B = [PP; ZP];                                        % all plankton

% Organic nitrogen
OM = reshape(v_in(fixedParams.OM_index), [nOM nz]);


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
att = (fixedParams.attSW + fixedParams.attP .* sum(PP)) .* fixedParams.zwidth';
att = 0.5 * att + [0 cumsum(att(1:nz-1))];
I = Isurf * exp(-att);


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
    QmaxVmax_ = params.Qmax_over_delQ .* params.Vmax_over_Qmin;
    aN = params.aN_over_Qmin .* N;
    mumax = params.Vmax_over_Qmin ./ (QmaxVmax_ ./ mu + 1); % maximum growth rate (1 / day) depends on nutrient uptake and photosynthetic rates    
    V = (gammaT .* mumax .* aN) ./ (aN + mumax); % growth rate
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

% Messy feeding
lambda_m = 1 - params.lambda_max;
beta_x_G = params.beta .* G;
ZP_x_lambda_m = lambda_m .* ZP;
OM_mess = ZP_x_lambda_m .* [sum(beta_x_G); sum(G - beta_x_G)]; % OM_mess=[DOM_mess; POM_mess]

% Mortality
beta_x_m_x_B = params.beta .* B_mortality;
OM_mort = [sum(beta_x_m_x_B); sum(B_mortality - beta_x_m_x_B)]; % OM_mort=[DOM_mort; POM_mort]

% Remineralisation
OM_remin = params.rOM .* OM;

SOM = OM_mort + OM_mess - OM_remin;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Remineralisation
SN = sum(OM_remin);

%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

v_diffuse = diffusion_1D([N(:), B', OM'], K, fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_diffuse = v_diffuse(:,2:nPP+2)';
OM_diffuse = v_diffuse(:,end-nOM+1:end)';

%~~~~~~~~
% Sinking
%~~~~~~~~

v_sink = sinking(OM', params.wk, fixedParams.zwidth);
OM_sink = v_sink';


%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - N_uptake_losses(:) + SN(:);
% Plankton
dBdt = B_diffuse + B_uptake + B_predation_gains - B_predation_losses - B_mortality;
% Organic matter
dOMdt = OM_diffuse + OM_sink + SOM;

% Separate phytoplankton and zooplankton
dPPdt = dBdt(fixedParams.phytoplankton,:);
dZPdt = dBdt(fixedParams.zooplankton,:);

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt(:)];


%% AUXILIARY OUTPUTS
if returnExtra
    % 1D
    extraOutput.PAR = I;
    
    % 2D (vectorised over cell size)
    if all(I == 0)
        extraOutput_2d.Qeq = repmat(params.Qmax, [1 nz]); % equilibrium nitrogen Qeq=Qmax when I=0
    else
        Q_ = (QmaxVmax_ .* aN) ./ (aN + params.Vmax_over_Qmin);
        extraOutput_2d.Qeq = (Q_ + mu) ./ (mu ./ params.Qmin + Q_ ./ params.Qmax); % equilibrium N quota
    end
    extraOutput_2d.PP_C = params.Q_C  .* (PP ./ extraOutput_2d.Qeq); % phytoplankton carbon biomass (mmol C / m^3)
end

end


%% Functions

function v = diffusion_1D(u,K,w,delz)
% Rate of change of u due to diffusion, assuming zero flux boundary conditions.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         K = diffusivities, size(K)=[nz-1 1]
%         w = depth layer widths, size(w)=[nz 1]
%         delz = distance between depth layer centers, size(delz)=[nz-1 1]
padZeros = zeros(1,size(u,2));
v = diff([padZeros; (K ./ delz) .* diff(u); padZeros]) ./ w;
end


function v = sinking(u,s,w)
% Rate of change of u due to sinking.
% Inputs: u = concentrations, size(u)=[nz nvar]
%         s = sinking speed, size(s)=[1 nvar]
%         w = depth layer widths, size(w)=[nz 1]
nvar = size(u,2);
v = ((-s) ./ w) .* diff([zeros(1,nvar); u]);
end





