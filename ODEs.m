function [dvdt, out1D, out2D, out3D] = ODEs(t, v_in, ... 
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
B = cat(1, ...
    reshape(v_in(fixedParams.PP_index), [nPP_size nz nPP_nut]), ... % autotrophs
    cat(3, ...
    reshape(v_in(fixedParams.ZP_index), [nZP_size nz nZP_nut]), ... % heterotrophs (include extra zeros for chl-a)
    zeros(nZP_size, nz, 1))); % all plankton

B_C = B(:,:,fixedParams.PP_C_index); % all planktonic carbon

% Organic matter
OM =reshape(v_in(fixedParams.OM_index), [nOM_type nz nOM_nut]);


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
att = (fixedParams.attSW + fixedParams.attP .* sum(B(phyto,:,fixedParams.PP_Chl_index))) .* fixedParams.zwidth';
% att = (fixedParams.attSW + fixedParams.attP .* sum(PP(:,:,fixedParams.PP_Chl_index))) .* fixedParams.zwidth';
att = 0.5 * att + [0 cumsum(att(1:nz-1))];
I = Isurf * exp(-att);


%% MODEL EQUATIONS

%~~~~~~~~~~~
% Physiology
%~~~~~~~~~~~

% nutrient quotas
Q = B ./ B_C;

% Nutrient limitation
gammaN = max(0, min(1, (Q(:,:,fixedParams.PP_N_index) - params.Qmin_QC) ./ params.delQ_QC));

% Uptake regulation
Qstat = 1 - gammaN .^ params.h;

% Temperature dependence
gammaT = exp(params.A .* (T - params.Tref));

% Background mortality
% B_C_mortality = params.m .* B_C; % linear mortality
% B_C_mortality = params.m .* B_C .^ 2; % non-linear mortality
mortality = params.m .* B;

%~~~~~~~~~~~
% Autotrophy
%~~~~~~~~~~~

V = zeros(nsize, nz, nPP_nut); % all uptake rates

% Nutrient uptake
V(:,:,fixedParams.PP_N_index) = ... 
    MichaelisMenton(params.Vmax_QC, params.kN, N) .* gammaT .* Qstat;
V(zoo,:,fixedParams.PP_N_index) = 0;

% Photosynthesis
zeroLight = all(I == 0);
if ~zeroLight
    psat = params.pmax .* gammaT .* gammaN; % light saturated photosynthetic rate
    aP_Q_I = (params.aP .* I) .* Q(:,:,fixedParams.PP_Chl_index);
    pc = psat .* (1 - exp(-aP_Q_I ./ psat )); % photosynthetic (carbon production) rate (1 / day)
    rho = params.theta .* pc ./ aP_Q_I;  % proportion of new nitrogen prodcution allocated to chlorophyll (mg Chl / mmol N)
    V(:,:,fixedParams.PP_Chl_index) = rho .* V(:,:,fixedParams.PP_N_index); % chlorophyll production rate (mg Chl / mmol C / day)
    V(:,:,fixedParams.PP_C_index) = max(0, pc - params.xi .* V(:,:,fixedParams.PP_N_index));
end

uptake = B_C .* V;
N_uptake_losses = sum(uptake(:,:,fixedParams.PP_N_index));


%~~~~~~~~~~~~~
% Heterotrophy
%~~~~~~~~~~~~~

% phi = exp(-log(fixedParams.delta ./ params.delta_opt) .^ 2 ./ (2 .* params.sigG .^ 2)); % phi could be moved outside of ODEs if delta_opt and sigG are not tuned
phi_BC = exp(-log(fixedParams.delta ./ params.delta_opt) .^ 2 ./ (2 .* params.sigG .^ 2)) .* ... 
    reshape(B_C, [1 size(B_C)]);
F = sum(phi_BC, 2);

phi_BC2 = phi_BC .^ 2;
% phi_BC2 = (phi .* reshape(B_C, [1 size(B_C)])) .^ 2;
Phi = phi_BC2 ./ sum(phi_BC2, 2); % prey preference

G = (reshape(gammaT, [1 size(gammaT)]) .* ...
    MichaelisMenton(params.Gmax, params.k_G, F) .* (1-exp(params.Lambda .* F))) .* Phi;  % grazing rate (1 / day)

predation_losses = reshape(Q, [1, nsize, nz, nPP_nut]) .* reshape(B_C(zoo,:), [nZP_size, 1, nz]) .* G;

lambda = zeros(nZP_size, 1, nz, nPP_nut);
lambda(:,:,:,fixedParams.ZP_C_index) = gammaN(zoo,:);
lambda(:,:,:,fixedParams.ZP_N_index) = Qstat(zoo,:);
lambda = params.lambda_max .* lambda;

predation_gains = lambda .* predation_losses;

mess = predation_losses(:,:,:,~fixedParams.PP_Chl_index) - ... 
    predation_gains(:,:,:,~fixedParams.PP_Chl_index);

if nZP_size > 1
    predation_losses = sum(predation_losses);  % sum over predators
end
predation_losses = reshape(predation_losses, [nsize, nz, nPP_nut]);
predation_gains = [zeros(nPP_size, nz, nPP_nut);  
    reshape(sum(predation_gains, 2), [nZP_size, nz, nPP_nut])];  % sum over prey


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of organic matter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Messy feeding
OM_mess = zeros(nOM_type, nz, nOM_nut);
beta_mess = reshape(params.beta, [1, nsize]) .* mess;
OM_mess(fixedParams.DOM_index,:,:) = sum(beta_mess, [1, 2]);
OM_mess(fixedParams.POM_index,:,:) = sum(mess - beta_mess, [1, 2]);

% Mortality
OM_mort = zeros(nOM_type, nz, nOM_nut);
beta_m_B = params.beta .* mortality(:,:,~fixedParams.PP_Chl_index);
OM_mort(fixedParams.DOM_index,:,:) = sum(beta_m_B);
OM_mort(fixedParams.POM_index,:,:) = sum(mortality(:,:,~fixedParams.PP_Chl_index) - beta_m_B);

% Remineralisation
OM_remin = params.rOM .* OM;

SOM = OM_mort + OM_mess - OM_remin; % (mmol / m^3 / day)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Remineralisation
SN = sum(OM_remin(:,:,fixedParams.OM_N_index)); % (mmol N / m^3 / day)

%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

B_C_t = B_C';
OM_ = reshape(permute(OM, [2, 1, 3]), [nz, nOM_type * nOM_nut]);

v_diffuse = diffusion_1D([N(:), B_C_t, OM_], K, fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_diffuse = Q .* v_diffuse(:,2:nsize+1)';
OM_diffuse = permute(reshape( ...
    v_diffuse(:,nsize+2:end), ...
    [nz, nOM_type, nOM_nut]), [2, 1, 3]);


%~~~~~~~~
% Sinking
%~~~~~~~~

wk = repmat(params.wk, [1 nOM_nut]);

v_sink = sinking([B_C_t, OM_], [params.wp, wk], fixedParams.zwidth);

B_sink = Q .* v_sink(:,1:nsize)';
OM_sink = permute(reshape(v_sink(:,nsize+1:end), ... 
    [nz, nOM_type, nOM_nut]), [2, 1, 3]);

%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - N_uptake_losses(:) + SN(:);

% Plankton
dBdt = B_sink + B_diffuse + uptake - predation_losses + predation_gains - mortality;
dPPdt = dBdt(phyto,:,:);
dZPdt = dBdt(zoo,:,~fixedParams.PP_Chl_index);

% Organic matter
dOMdt = OM_diffuse + OM_sink + SOM;

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt(:)];



%% AUXILIARY OUTPUTS

if (islogical(returnExtra) && returnExtra) || ... 
        (~islogical(returnExtra) && ~any(strcmp(returnExtra, 'none')))
    
    %~~~~~~~~~~~
    % 1D (depth)
    %~~~~~~~~~~~
    out.PAR = I; % PAR-at-depth depends upon plankton concentrations
    out.gammaT = gammaT; % Temperature dependence
    
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % 2D (depth & cell size, including zooplankton)
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out.cellDensity = B_C ./ params.Q_C;
    out.biovolume = 1e-18 * fixedParams.sizeAll .* out.cellDensity;
    out.gammaN = gammaN;
    out.Qstat = Qstat;
    
    if ~zeroLight
        out.psat = psat;
        out.pc = pc;
        out.rho = rho;
    else
        out.psat = zeros_size_nz;
        out.pc = zeros_size_nz;
        out.rho = zeros_size_nz;
    end
    
    % 3D (depth & size & nutrient)
    out.Q = Q;
    out.V = V;
    out.predation_losses = predation_losses;
    out.predation_gains = predation_gains;
    
    if ~islogical(returnExtra) && ~any(strcmp(returnExtra,'all'))
        f = fieldnames(out);
        f = f(~ismember(f, returnExtra));
        out = rmfield(out, f);
        if isempty(out), clear out; end
    end
    
    [out1D, out2D, out3D] = groupByDimension(out);
    
else
    out1D = []; out2D = []; out3D = [];
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

function [out1, out2, out3] = groupByDimension(v)
% Organises extra output structs
fields = fieldnames(v);
out1 = struct();
out2 = struct();
out3 = struct();
for i = 1:length(fields)
    x = v.(fields{i});
    if isvector(x)
        out1.(fields{i}) = x;
    end
    if ~isvector(x) && ismatrix(x)
        out2.(fields{i}) = x;
    end
    if size(x, 1) > 1 && ndims(x) == 3
        out3.(fields{i}) = x;
    end
end
end

