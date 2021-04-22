function [dvdt, out] = ODEs(t, v_in, parameterList, forc, timeStep, returnExtra)

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
out.I = Isurf * exp(-att);


%% MODEL EQUATIONS

%~~~~~~~~~~~
% Physiology
%~~~~~~~~~~~

% nutrient quotas
out.Q = B ./ B_C;

% Nutrient limitation
out.gammaN = max(0, min(1, (out.Q(:,:,fixedParams.PP_N_index) - params.Qmin_QC) ./ params.delQ_QC));

% Uptake regulation
out.Qstat = 1 - out.gammaN .^ params.h;

% Temperature dependence
out.gammaT = exp(params.A .* (T - params.Tref));

% Background mortality
out.mortality = params.m .* B;

%~~~~~~~~~~~
% Autotrophy
%~~~~~~~~~~~

out.V = zeros(nsize, nz, nPP_nut); % all uptake rates

% Nutrient uptake
out.V(:,:,fixedParams.PP_N_index) = ... 
    MichaelisMenton(params.Vmax_QC, params.kN, N) .* out.gammaT .* out.Qstat;
out.V(zoo,:,fixedParams.PP_N_index) = 0;

% Photosynthesis
zeroLight = all(out.I == 0);
if ~zeroLight
    out.psat = params.pmax .* out.gammaT .* out.gammaN; % light saturated photosynthetic rate
    aP_Q_I = (params.aP .* out.I) .* out.Q(:,:,fixedParams.PP_Chl_index);
    out.pc = out.psat .* (1 - exp(-aP_Q_I ./ out.psat )); % photosynthetic (carbon production) rate (1 / day)
    out.rho = params.theta .* out.pc ./ aP_Q_I;  % proportion of new nitrogen prodcution allocated to chlorophyll (mg Chl / mmol N)
    out.V(:,:,fixedParams.PP_Chl_index) = out.rho .* out.V(:,:,fixedParams.PP_N_index); % chlorophyll production rate (mg Chl / mmol C / day)
    out.V(:,:,fixedParams.PP_C_index) = max(0, out.pc - params.xi .* out.V(:,:,fixedParams.PP_N_index));
end

out.uptake = B_C .* out.V;
out.N_uptake_losses = sum(out.uptake(:,:,fixedParams.PP_N_index));


%~~~~~~~~~~~~~
% Heterotrophy
%~~~~~~~~~~~~~

phi_BC = params.phi .* reshape(B_C, [1 size(B_C)]);
F = sum(phi_BC, 2);
phi_BC2 = phi_BC .^ 2;
Phi = phi_BC2 ./ sum(phi_BC2, 2); % prey preference

out.G = (reshape(out.gammaT, [1 size(out.gammaT)]) .* ...
    MichaelisMenton(params.Gmax, params.k_G, F) .* (1-exp(params.Lambda .* F))) .* Phi;  % grazing rate (1 / day)

out.predation_losses_all = reshape(out.Q, [1, nsize, nz, nPP_nut]) .* reshape(B_C(zoo,:), [nZP_size, 1, nz]) .* out.G;

out.lambda = zeros(nZP_size, 1, nz, nPP_nut);
out.lambda(:,:,:,fixedParams.ZP_C_index) = out.gammaN(zoo,:);
out.lambda(:,:,:,fixedParams.ZP_N_index) = out.Qstat(zoo,:);
out.lambda = params.lambda_max .* out.lambda;

out.predation_gains_all = out.lambda .* out.predation_losses_all;

mess = out.predation_losses_all(:,:,:,~fixedParams.PP_Chl_index) - ... 
    out.predation_gains_all(:,:,:,~fixedParams.PP_Chl_index);

if nZP_size > 1
    out.predation_losses = sum(out.predation_losses_all);  % sum over predators
end
out.predation_losses = reshape(out.predation_losses, [nsize, nz, nPP_nut]);
out.predation_gains = [zeros(nPP_size, nz, nPP_nut);  
    reshape(sum(out.predation_gains_all, 2), [nZP_size, nz, nPP_nut])];  % sum over prey


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources and sinks of organic matter
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Messy feeding
out.OM_mess = zeros(nOM_type, nz, nOM_nut);
beta_mess = reshape(params.beta, [1, nsize]) .* mess;
out.OM_mess(fixedParams.DOM_index,:,:) = sum(beta_mess, [1, 2]);
out.OM_mess(fixedParams.POM_index,:,:) = sum(mess - beta_mess, [1, 2]);

% Mortality
out.OM_mort = zeros(nOM_type, nz, nOM_nut);
beta_m_B = params.beta .* out.mortality(:,:,~fixedParams.PP_Chl_index);
out.OM_mort(fixedParams.DOM_index,:,:) = sum(beta_m_B);
out.OM_mort(fixedParams.POM_index,:,:) = sum(out.mortality(:,:,~fixedParams.PP_Chl_index) - beta_m_B);

% Remineralisation
out.OM_remin = params.rOM .* OM;

SOM = out.OM_mort + out.OM_mess - out.OM_remin; % (mmol / m^3 / day)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Sources of inorganic nutrients
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Remineralisation
SN = sum(out.OM_remin(:,:,fixedParams.OM_N_index)); % (mmol N / m^3 / day)

%~~~~~~~~~~
% Diffusion
%~~~~~~~~~~

B_C_t = B_C';
OM_ = reshape(permute(OM, [2, 1, 3]), [nz, nOM_type * nOM_nut]);

v_diffuse = diffusion_1D([N(:), B_C_t, OM_], K, fixedParams.zwidth, fixedParams.delz);

N_diffuse = v_diffuse(:,1);
B_diffuse = out.Q .* v_diffuse(:,2:nsize+1)';
OM_diffuse = permute(reshape( ...
    v_diffuse(:,nsize+2:end), ...
    [nz, nOM_type, nOM_nut]), [2, 1, 3]);


%~~~~~~~~
% Sinking
%~~~~~~~~

wk = repmat(params.wk, [1 nOM_nut]);

v_sink = sinking([B_C_t, OM_], [params.wp, wk], fixedParams.zwidth);

B_sink = out.Q .* v_sink(:,1:nsize)';
OM_sink = permute(reshape(v_sink(:,nsize+1:end), ... 
    [nz, nOM_type, nOM_nut]), [2, 1, 3]);

%~~~~~
% ODEs
%~~~~~

% Inorganic nutrients
dNdt = N_diffuse - out.N_uptake_losses(:) + SN(:);

% Plankton
dBdt = B_sink + B_diffuse + out.uptake - out.predation_losses + out.predation_gains - out.mortality;
dPPdt = dBdt(phyto,:,:);
dZPdt = dBdt(zoo,:,~fixedParams.PP_Chl_index);

% Organic matter
dOMdt = OM_diffuse + OM_sink + SOM;

dvdt = [dNdt; dPPdt(:); dZPdt(:); dOMdt(:)];



%% AUXILIARY OUTPUTS

if (islogical(returnExtra) && returnExtra) || ... 
        (~islogical(returnExtra) && ~any(strcmp(returnExtra, 'none')))
    
    out.cellDensity = B_C ./ params.Q_C;
    out.biovolume = 1e-18 * fixedParams.sizeAll .* out.cellDensity;

    % Extra output variables retained by default when return = true or 'all'.
    keepVars = {'I', 'Q', 'V', 'G', 'lambda', 'cellDensity', 'biovolume'};
    % Any term can be included in keepVars, but it's useful to be sparing
    % with memory by only on returning some terms then deriving extra output
    % outside the ODEs.m script.
    
    if ~islogical(returnExtra) && ~all(strcmp(returnExtra, 'all'))
        % if extra output variables have been specified explicitly...
        keepVars = returnExtra;
    end
    
    fields = fieldnames(out);
    out = rmfield(out, fields(~ismember(fields, keepVars)));
    
else
    out = struct();
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

% function [out1, out2, out3, out4] = groupByDimension(v)
% % Organises extra output structs
% fields = fieldnames(v);
% out1 = struct();
% out2 = struct();
% out3 = struct();
% out4 = struct();
% for i = 1:length(fields)
%     x = v.(fields{i});
%     if isvector(x)
%         out1.(fields{i}) = x;
%     end
%     if ~isvector(x) && ismatrix(x)
%         out2.(fields{i}) = x;
%     end
% %     if size(x, 1) > 1 && ndims(x) == 3
%     if ndims(x) == 3
%         out3.(fields{i}) = x;
%     end
% %     if size(x, 1) > 1 && ndims(x) == 4
%     if ndims(x) == 4
%         out4.(fields{i}) = x;
%     end
% end
% end

