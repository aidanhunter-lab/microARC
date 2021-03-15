function [dvdt, extraOutput, extraOutput_2d] = ODEs(t, v_in, ... 
    parameterList, forc, timeStep, returnExtra)

fixedParams = parameterList.FixedParams;
params = parameterList.Params;

%% MODEL DIMENSIONS

nz = fixedParams.nz;    % number of depth layers
nPP_size = fixedParams.nPP_size;  % number of phytoplankton size classes
nPP_nut = fixedParams.nPP_nut;  % number of phytoplankton nutrient classes
nZP_size = fixedParams.nZP_size;  % number of zooplankton size classes
nZP_nut = fixedParams.nZP_nut;  % number of zooplankton nutrient classes
nOM_type = fixedParams.nOM_type;  % number of organic matter types
nOM_nut = fixedParams.nOM_nut;  % number of organic nutrient classes
phyto = fixedParams.phytoplankton;
zoo = fixedParams.zooplankton;
nsize = nPP_size + nZP_size;

%% INITIAL CONDITIONS

% Inorganic nitrogen
N = v_in(fixedParams.IN_index)';

% Plankton
PP = reshape(v_in(fixedParams.PP_index), [nPP_size nz nPP_nut]); % phytoplankton (all nutrients)
% P_C = PP(:,:,fixedParams.PP_C_index);
ZP = reshape(v_in(fixedParams.ZP_index), [nZP_size nz nZP_nut]);
% Z_C = ZP(:,:,fixedParams.ZP_C_index); % zooplankton (carbon)

B = cat(1, PP, ...
    cat(3, ZP, zeros(1, nz, 1))); % all plankton

B_C = B(:,:,fixedParams.PP_C_index); % all planktonic carbon


% Organic matter
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

% N and Chl quotas relative to C
Q_N = B(:,:,fixedParams.PP_N_index) ./ B_C;
Q_Chl = B(:,:,fixedParams.PP_Chl_index) ./ B_C;

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

% zeros_size_nz = zeros(nPP_size, nz);
zeros_size_nz = zeros(nsize, nz);

% Nutrient uptake
V_N = MichaelisMenton(params.Vmax_QC, params.kN, N) .* gammaT .* Qstat;
V_N(zoo,:) = 0;

% V_N = params.Vmax_QC ./ (1 + params.Vmax_QC ./ (params.aN_QC .* N)) .* gammaT .* Qstat; % nitrogen uptake rate (mmol N / mmol C / day)
N_uptake = V_N .* B_C; % mmol N / m^3 / day
N_uptake_losses = sum(N_uptake);

% Photosynthesis
zeroLight = all(I == 0);
if zeroLight
    V_Chl = zeros_size_nz;
    V_C = zeros_size_nz;
%     V_C = zeros(nPP_size+1, nz);
else
    psat = params.pmax .* gammaT .* gammaN; % light saturated photosynthetic rate
    aP_Q_I = (params.aP .* I) .* Q_Chl;
    pc = psat .* (1 - exp(-aP_Q_I ./ psat )); % photosynthetic (carbon production) rate (1 / day)
    rho = params.theta .* pc ./ aP_Q_I;  % proportion of new nitrogen prodcution allocated to chlorophyll (mg Chl / mmol N)
    V_Chl = rho .* V_N; % chlorophyll production rate (mg Chl / mmol C / day)
    V_C = max(0, pc - params.xi .* V_N);
%     V_C = [max(0, pc - params.xi .* V_N); zeros(1, nz)];
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

predation_losses_C = G .* B_C; % mmol C / m^3 / day
predation_gains_C = [zeros(nPP_size, nz); ... 
    params.lambda_max .* sum(predation_losses_C)];
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of organic matter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Messy feeding
lambda_predLoss_C = (1 - params.lambda_max) .* predation_losses_C;
beta_lambda_predLoss_C = params.beta .* lambda_predLoss_C;
OM_mess_C = [sum(beta_lambda_predLoss_C); ... 
    sum(lambda_predLoss_C-beta_lambda_predLoss_C)]; % OM_mess=[DOM_mess; POM_mess] (mmol C / m^3 / day)

lambda_predLoss_N = Q_N .* lambda_predLoss_C;
beta_lambda_predLoss_N = params.beta .* lambda_predLoss_N;
OM_mess_N = [sum(beta_lambda_predLoss_N); ... 
    sum(lambda_predLoss_N-beta_lambda_predLoss_N)]; % OM_mess=[DOM_mess; POM_mess] (mmol N / m^3 / day)

% Mortality
beta_m_B = params.beta .* B_C_mortality;
OM_mort_C = [sum(beta_m_B); sum(B_C_mortality - beta_m_B)]; % OM_mort=[DOM_mort; POM_mort] (mmol C / m^3 / day)

B_N_mortality = Q_N .* B_C_mortality;
beta_m_B = params.beta .* B_N_mortality;
OM_mort_N = [sum(beta_m_B); sum(B_N_mortality - beta_m_B)]; % OM_mort=[DOM_mort; POM_mort] (mmol N / m^3 / day)

% Remineralisation
OM_remin = params.rOM .* OM;

SOM_C = OM_mort_C + OM_mess_C - OM_remin(:,:,fixedParams.OM_C_index); % (mmol C / m^3 / day)
SOM_N = OM_mort_N + OM_mess_N - OM_remin(:,:,fixedParams.OM_N_index); % (mmol N / m^3 / day)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Remineralisation
SN = sum(OM_remin(:,:,fixedParams.OM_N_index)); % (mmol N / m^3 / day)

%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

B_C_t = B_C';
OM_C_t = OM_C';

v_diffuse = diffusion_1D([N(:), B_C_t, OM_C_t], K, fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_C_diffuse = v_diffuse(:,2:nPP_size+nZP_size+1)';
OM_C_diffuse = v_diffuse(:,end-nOM_type+1:end)';

%~~~~~~~~
% Sinking
%~~~~~~~~

v_sink = sinking([B_C_t, OM_C_t], [params.wp, params.wk], fixedParams.zwidth);
% v_sink = sinking([B_C_t(:,phyto), OM_C_t], [params.wp, params.wk], fixedParams.zwidth);

B_C_sink = v_sink(:,1:nPP_size+nZP_size)';
% B_C_sink = [v_sink(:,1:nPP_size), zeros(nz, 1)]';
OM_C_sink = v_sink(:,nPP_size+nZP_size+1:end)';

%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - N_uptake_losses(:) + SN(:);

% Plankton
fluxC_ = B_C_sink + B_C_diffuse - predation_losses_C - B_C_mortality; % mmol C / m^3 / day
fluxC = fluxC_ + predation_gains_C;
fluxN = N_uptake + Q_N .* fluxC; % mmol N / m^3 / day
fluxChl = V_Chl .* B_C + Q_Chl .* fluxC_; % mg Chl / m^3 / day
fluxC = V_C .* B_C + fluxC;
dPPdt = cat(3, fluxC(phyto,:), fluxN(phyto,:), fluxChl(phyto,:));
dZPdt = cat(3, fluxC(zoo,:), fluxN(zoo,:));

% Organic matter
fluxC_ = OM_C_diffuse + OM_C_sink;
fluxN = OM(:,:,fixedParams.OM_N_index) ./ OM_C .* fluxC_ + SOM_N;
fluxC = fluxC_ + SOM_C;
dOMdt = cat(3, fluxC, fluxN);

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt(:)];



%% AUXILIARY OUTPUTS

if (islogical(returnExtra) && returnExtra) || ... 
        (~islogical(returnExtra) && ~any(strcmp(returnExtra, 'none')))
    
    cellDensity = B_C ./ params.Q_C;
    biovolume = 1e-18 * fixedParams.sizeAll .* cellDensity;
    
    %~~~~~~~~~~~
    % 1D (depth)
    %~~~~~~~~~~~
    extraOutput.PAR = I; % PAR-at-depth depends upon plankton concentrations
    extraOutput.gammaT = gammaT; % Temperature dependence
    extraOutput.F = F; % Total prey carbon
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 2D (depth & cell size, including zooplankton)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    extraOutput_2d.cellDensity = cellDensity;
    extraOutput_2d.biovolume = biovolume;
    extraOutput_2d.Q_N = Q_N;
    extraOutput_2d.Q_Chl = Q_Chl;
    extraOutput_2d.gammaN = gammaN;
    extraOutput_2d.Qstat = Qstat;
    extraOutput_2d.V_N = V_N;
    
    if ~zeroLight
        extraOutput_2d.V_Chl = V_Chl;
        extraOutput_2d.V_C = V_C;
        extraOutput_2d.psat = psat;
        extraOutput_2d.pc = pc;
        extraOutput_2d.rho = rho;
    else
        extraOutput_2d.V_Chl = zeros_size_nz;
        extraOutput_2d.V_C = zeros_size_nz;
        extraOutput_2d.psat = zeros_size_nz;
        extraOutput_2d.pc = zeros_size_nz;
        extraOutput_2d.rho = zeros_size_nz;
    end
    extraOutput_2d.Phi = Phi;
    extraOutput_2d.G = G;
    extraOutput_2d.predation_losses_C = predation_losses_C;
    extraOutput_2d.predation_gains_C = predation_gains_C;
    
%     fields = fieldnames(extraOutput_2d);
%     for i = 1:length(fields)
%         if size(extraOutput_2d.(fields{i}), 1) == nPP_size
%             extraOutput_2d.(fields{i}) = [extraOutput_2d.(fields{i}); nan(1, nz)];
%         end
%     end
    
    
    
    if ~islogical(returnExtra) && ~any(strcmp(returnExtra,'all'))
        f1 = fieldnames(extraOutput);
        f1 = f1(~ismember(f1, returnExtra));
        f2 = fieldnames(extraOutput_2d);
        f2 = f2(~ismember(f2, returnExtra));
        extraOutput = rmfield(extraOutput, f1);
        extraOutput_2d = rmfield(extraOutput_2d, f2);
        if isempty(extraOutput), clear extraOutput; end
        if isempty(extraOutput_2d), clear extraOutput_2d; end
    end
    
    
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

