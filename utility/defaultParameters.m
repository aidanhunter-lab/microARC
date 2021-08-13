function [FixedParams, Params, Bounds] = defaultParameters(bioModel)
% Set default parameter values.
% Some variables left empty will be assigned values in initialiseParameters.m

%% Fixed parameters

FixedParams.bioModel = bioModel;

%~~~~~~~~
% Timings
%~~~~~~~~
FixedParams.nt    = []; % number of modelled time steps (determined from forcing data)
FixedParams.years = []; % number of modelled years (determined from forcing and fitting data)

%~~~~~~~~~~~~~
% Depth layers
%~~~~~~~~~~~~~
% Modelled layer width increases with depth => model has greater resolution
% near surface waters. Define 2 depth zones: from the surface to depth
% Hmain; and from depth Hmain to Htot. The 1st of these zones should be set
% where the interesting plankton dynamics occurs -- where there's
% sufficient light. Within this zone, the depth layer widths increase with
% depth. The 2nd zone is included to model variables at depths below the 
% nitrocline, and should be set with reference to nutrient data. This zone
% is simply a single extra wide modelled layer.
FixedParams.nz     = 9;                                                    % number of modelled depth layers
FixedParams.Hmain = 100;
FixedParams.Htot  = 300;                                                   % total modelled depth
FixedParams.dzmin  = 5;                                                    % minimum layer width (dzmin <= Htot/(nz-1))
FixedParams.dzmax  = 2 * FixedParams.Hmain / (FixedParams.nz - 1) - ... 
    FixedParams.dzmin;                                                     % maximum layer width
FixedParams.zwidth = linspace(FixedParams.dzmin, FixedParams.dzmax, ... 
    FixedParams.nz - 1)';                                                  % widths of depth layers
FixedParams.zwidth = [FixedParams.zwidth; ... 
    FixedParams.Htot - FixedParams.Hmain];
FixedParams.zw     = [0; -cumsum(FixedParams.zwidth)];                     % depth of layer edges
FixedParams.z      = [];                                                   % midpoints of depth layers
FixedParams.delz   = [];                                                   % distance between centres of adjacent depth layers

%~~~~~~~~~~~~~~~~
% State variables
%~~~~~~~~~~~~~~~~

% Inorganic nutrient
FixedParams.IN_nut = {'NO3'}; % nutrient types - only nitrate is modelled
FixedParams.nIN    = [];

% Plankton
FixedParams.PP_nut  = {'C','N','Chl'};                                % phytoplankton cell contents - carbon, nitrogen and chlorophyll
FixedParams.ZP_nut  = {'C','N'};                                      % zooplankton cell contents - carbon, nitrogen
FixedParams.nPP_nut = [];
FixedParams.nZP_nut = [];

FixedParams.nPP_size = 9; % number of phytoplankton size classes
FixedParams.PPdia_intervals = [];
FixedParams.PPdia = 2 .^ (0:FixedParams.nPP_size - 1); % cell diameters -- each class is double the diameter of the previous (evenly spaced on log scale)
FixedParams.PPdia = FixedParams.PPdia(:);
logPPdia = log10(FixedParams.PPdia);
d = diff(logPPdia(1:2));
FixedParams.PPdia_intervals = 10 .^ [logPPdia - 0.5 .* d; ... 
    logPPdia(end) + 0.5 .* d];   % edges of size class intervals
FixedParams.PPsize = d2vol(FixedParams.PPdia); % cell volumes [mu m^3]

% Use the same size interval grid for heterotrophs
FixedParams.nZP_size = FixedParams.nPP_size;
FixedParams.ZPdia_intervals = FixedParams.PPdia_intervals;
FixedParams.ZPdia = FixedParams.PPdia;
FixedParams.ZPsize = FixedParams.PPsize;

FixedParams.nPP           = [];                                       % number of phytoplankton classes
FixedParams.nZP           = [];                                       % number of heterotroph classes

FixedParams.diaAll = [FixedParams.PPdia; FixedParams.ZPdia];
FixedParams.sizeAll = [FixedParams.PPsize; FixedParams.ZPsize];

% delta_pred = reshape(FixedParams.ZPsize, [FixedParams.nZP_size, 1]);
% delta_prey = [reshape(FixedParams.PPsize, [1 FixedParams.nPP_size]), ... 
%     reshape(FixedParams.ZPsize, [1 FixedParams.nZP_size])];
% FixedParams.delta = delta_pred ./ delta_prey; % predator:prey size ratios
FixedParams.delta = []; % predator:prey size ratios

FixedParams.diatoms       = [];                                       % assume large phytoplankton are diatoms - only needed to split SINMOD output over classes during state variable initialisation
FixedParams.phytoplankton = [];                                       % index phytoplankton
FixedParams.zooplankton   = [];                                       % index zooplankton

% Organic matter
FixedParams.OM_nut   = {'C', 'N'};
FixedParams.nOM_nut  = [];
FixedParams.OM_type  = {'DOM', 'POM'};
FixedParams.nOM_type = [];
FixedParams.nOM      = [];

% All variables
FixedParams.nVar       = [];                                          % number of state variables per depth and trajectory
FixedParams.nEquations = [];                                          % total number of ODEs

%~~~~~~~~~~~~~~~~~~~~~~~~~~~
% All other fixed quantities
%~~~~~~~~~~~~~~~~~~~~~~~~~~~
FixedParams.m_min = 0.01;
FixedParams.attSW = 0.04;                                  % light attenuation by sea water (1 / m)
FixedParams.attP = 0.0149;                                 % light attenuation by chlorophyll (m^2 / mg Chl), Krause-Jensen & Sand-Jensen, Limnol. Oceanogr. 43(3):396-407(1998)
FixedParams.minSinkSpeed_POM = 1;                          % minimum sinking speed (m / day) of POM (this is used to constrain POM sink speed if it's assumed to be depth dependent -- slower speeds near the surface fit the data better and may be justified if choppyness near the surface increases POM hang-time)
FixedParams.maxSinkSpeed_POM = 10;                         % maximum sinking speed (m / day) of POM (this determines max integration time step)
FixedParams.depthDependentPOMsinkSpeed = false;            % can set to true to model POM sink speed as depth-dependent, with lower speeds near the surface
maxDepthLayerClearRate = 0.1;                              % max rate at which depth layer can by cleared of plankton via sinking -- max sink speed of plankton is specified like this for model-stability purposes
FixedParams.maxSinkSpeed_P = maxDepthLayerClearRate .* ...
    min(FixedParams.zwidth);                               % max sinking speed of plankton (when this is too great the model becomes unstable because cells clear from depth layers too quickly)
FixedParams.POM_is_lost  = true;                           % is POM lost from the system by sinking below bottom modelled depth layer
FixedParams.NclineDepth = 110;                             % depth of N-cline -- in-situ data sampled below NclineDepth should not vary much with depth, check plots of data to set this parameter.

%% Variable parameters

% List names of all parameters that MAY be numerically optimised
Params.scalarParams = {
    'Tref'
    'A'
    'h'
    'm2'
    'aP'
    'theta'
    'xi'
    'k_G'
    'delta_opt'
    'sigG'
    'Lambda'
    'lambda_max'
    'rDON'
    'rDOC'
    'rPON'
    'rPOC'
    };
Params.vectorParams = {
    % These are scalar parameters that are combined in functions to create vectors
    'Q_C_a'
    'Q_C_b'
    'Qmin_QC_a'
    'Qmin_QC_b'
    'Qmax_delQ_a'
    'Qmax_delQ_b'
    'Vmax_QC_a'
    'Vmax_QC_b'
    'aN_QC_a'
    'aN_QC_b'
    'pmax_a'
    'pmax_b'
    'Gmax_a'
    'Gmax_b'
    'm_a'
    'm_b'
    'K_m_coef'
    'beta1'
    'beta2'
    'beta3'
    'wp_a'
    'wp_b'
    'wDOM1'
    'wPOM1'
    };


%~~~~~~~~~
% Plankton
%~~~~~~~~~

% Size-dependent
% (mostly power functions y=aV^b)

% carbon quota [mmol C / cell], Maranon et al. (2013)
Params.Q_C_func = @(a,b,V) powerFunction(a,b,V);
Params.Q_C_a = (1/12) * 1e-9 * 10^-0.69;
Params.Q_C_b = 0.88;
Params.Q_C = [];

% min and max cellular nitrogen quota [mmol N / cell], scaled by C quota [mmol N / mmol C], Maranon et al. (2013)
Params.Qmin_QC_func = @(a,b,V) powerFunction(a,b,V);
Params.Qmin_QC_a = (1/14) * 1e-9 * 10^-1.47 / Params.Q_C_a;
% Params.Qmin_QC_b = 0.84 - Params.Q_C_b;
Params.Qmin_QC_b = 0; % minimum N:C ratio expected to have same size dependence as C-quota => set Qmin_QC_b = 0
Params.Qmin_QC = [];

% max quota

Params.Qmax_QC_func = @(Qmin_QC, Qmax_delQ) Qmin_QC ./ (1 - 1 ./ Qmax_delQ);
Params.Qmax_QC = [];

Params.delQ_QC_func = @(Qmin_QC, Qmax_QC) Qmax_QC - Qmin_QC;
Params.delQ_QC = [];

Params.Qmax_delQ_func = @(a,b,V) 1 ./ (1 - (a .* V .^ b));
% Params.Qmax_delQ_func = @(a,b,V) 1 ./ (1 - (a .* V .^ -exp(b))); % estimate b on negative log scale
Params.Qmax_delQ_a = 10^(-1.47+1.26);
Params.Qmax_delQ_b = 0.84 - 0.93;
% Params.Qmax_delQ_b = log(-Params.Qmax_delQ_b); % estimate on log negative scale
Params.Qmax_delQ = [];

% maximum nitrogen uptake rate Vmax [mmol N/cell/day] scaled by C quota, Vmax_over_QC [mmol N / mmol C / day], Maranon et al. (2013)
Params.Vmax_QC_func = @(a,b,V) powerFunction(a,b,V);
Params.Vmax_QC_a = 24 / 14 * 1e-9 * 10^-3 / Params.Q_C_a;
Params.Vmax_QC_b = 0.97 - Params.Q_C_b;
Params.Vmax_QC = [];

% cellular affinity for nitrogen scaled by QC [m^3 / mmol C / day], derived using half saturation from Litchmann et al. (2007)
Params.aN_QC_func = @(a,b,V) powerFunction(a,b,V);
% Params.aN_QC_func = @(a,b,V) powerFunction(a,-exp(b),V); % estimate b on negative log scale
Params.aN_QC_a = Params.Vmax_QC_a / 10^-0.77;
Params.aN_QC_b = Params.Vmax_QC_b -0.27;
% Params.aN_QC_b = log(-Params.aN_QC_b); % estimate on log negative scale
Params.aN_QC = [];

% half saturation
Params.kN_func = @(Vmax_QC, aN_QC) Vmax_QC ./ aN_QC;
Params.kN = [];

% maximum photosynthetic rate [1/day]
Params.pmax_func = @(a,b,V) powerFunction(a,b,V);
% Params.pmax_func = @(a,b,V) powerFunction(a,-exp(b),V); % b estimated on negative scale
Params.pmax_a = 2.5;
Params.pmax_b = -0.15;
% Params.pmax_b = log(-Params.pmax_b); % estimate on negative log scale
Params.pmax = [];

% maximum grazing rate
Params.Gmax_func = @(a,b,V) powerFunction(a,b,V);
% Params.Gmax_func = @(a,b,V) powerFunction(a,-exp(b),V); % b estimated on negative scale
Params.Gmax_a = 21;
Params.Gmax_b = -0.16;
% Params.Gmax_b = log(-Params.Gmax_b); % estimate on negative log scale
Params.Gmax = [];

% % half-saturation prey concentration (mmol N / m^3) for grazing uptake
% Params.k_G_a = 0.5;
% Params.k_G_b = 0.18;
% Params.k_G = [];

% prey size preferences
Params.phi_func = @(delta_opt, sigG, delta) ... 
    exp(-log(delta ./ delta_opt) .^ 2 ./ (2 .* sigG .^ 2));
Params.phi = [];

% background mortality -- linear
Params.m_func = @(a,b,V) powerFunction(a,b,V);
% Params.m_func = @(a,b,V) powerFunction(a,-exp(b),V);
Params.m_a = 0.05; % mortality for cell volume = 1 mu m ^ 3
Params.m_b = 0; % mortality size-exponent
% Params.m_b = -0.1; % mortality size-exponent
% Params.m_b = log(-Params.m_b); % estimate on negative log scale
Params.m = [];

% Hyperbolic term in mortality expression may be used to prevent extremely
% low abundances/extinction during prolonged poor growth conditions.
% Setting K_m_coef=0 removes this effect.
Params.K_m_func = @(K_m_coef, Q_C) K_m_coef .* Q_C; % I'm not convinced that multiplying Q_C is a good approach... but just now I cant think of a better and principled method...
% Params.K_m_coef = 1; % Params.K_m_coef = 1 is equivalent to mortality 'half-saturation' at concentration of 1 cell / m^3
Params.K_m_coef = 0;
Params.K_m = [];

% sinking plankton (parameter sensitivity assessment suggests omitting
% plankton sinking => set to zero)
Params.wp_func = @(a,b,V) powerFunction(a,b,V);
% Params.wp_a = 1e-5; % intercept = 0 => no sinking at any size
Params.wp_a = 0; % intercept = 0 => no sinking at any size
% Params.wp_b = 0.33;
Params.wp_b = 0;
Params.wp = [];

% partitioning of dead matter into DOM and POM,  initial values from Ward et al. (2016)
Params.beta_func = @(b1,b2,b3,V) b1 ./ (1 + exp(log10(V) - b3)) + ...
    b1 .* b2 ./ (1 + exp((b3 - log10(V)))); % % flexible 3-parameter double logistic function of log10(V)
Params.beta1 = 0.9;
Params.beta2 = 0.2/0.9;
Params.beta3 = 2.0;
Params.beta = [];

% Size-independent
% Params.m2 = 1e-5;           % quadratic, density-dependent, term in background mortality (1/day)
Params.m2 = 0;              % param sensitivity analysis suggests setting quadratic mortality term to zero...
Params.aP = 3.83e-7;        % initial slope of photosynthesis-irradiance curve [(mmol C m^2) / (mg Chl mu Ein)]
Params.theta = 4.2;         % maximum chl:N ratio [mg Chl / mmol N], Geider et al., 1998
Params.xi = 2.33;           % cost of photosynthesis (mmol C / mmol N)
Params.Tref = 20;           % reference temperature (degrees C)
Params.A = 0.05;            % temperature dependence (unitless)
Params.h = 10;              % curvature on quota uptake limitation
% Params.m = 0.05;            % linear plankton mortality (1/day)
Params.k_G = 5;             % half-saturation prey concentration (mmol C / m^3) for grazing uptake
% Params.Gmax = 5;            % maximum grazing rate
Params.delta_opt = 10;      % optimal predator:prey size ratio maximising feeding fluxes
if strcmp(bioModel, 'singlePredatorClass')
    % adjust size ratios for single predator class model -- all prey grazed equally
    FixedParams.delta = repmat(Params.delta_opt, size(FixedParams.delta));
end
Params.sigG = 2;            % variability of predator:prey size ratios
Params.Lambda = -1;         % prey refuge parameter (unitless)
Params.lambda_max = 0.7;    % maximum prey assimilation efficiency

%~~~~~~~~~~~~~~~
% Organic matter
%~~~~~~~~~~~~~~~
Params.rOM_func = @(rDOC, rDON, rPOC, rPON, DOM_ind, POM_ind, C_ind, N_ind) ...
    reshape( ...
    rDOC .* (DOM_ind' & C_ind) + ...
    rDON .* (DOM_ind' & N_ind) + ...
    rPOC .* (POM_ind' & C_ind) + ...
    rPON .* (POM_ind' & N_ind), ...
    [length(DOM_ind), 1 , length(C_ind)]);
Params.rOM = []; % remineralisation rates (1 / day)

% Assume that remineralisation rates are identical for all nutrients within
% each OM type, i.e., choose values for N, then set values for C identically
Params.rDOC_func = @(rDON) rDON;
Params.rDON = 0.02;         % dissolved organic nitrogen remineralisation rate (1/day)
% Params.rDOC = 0.02;         % dissolved organic carbon remineralisation rate (1/day)
Params.rDOC = [];           % dissolved organic carbon remineralisation rate (1/day)
Params.rPOC_func = @(rPON) rPON;
Params.rPON = 0.04;         % particulate organic nitrogen remineralisation rate (1/day)
Params.rPOC = [];           % particulate organic nitrogen remineralisation rate (1/day)

% Sinking: parameters wDOM1 and wPOM1 are scalars that, together with
% functions wDOM_func and wPOM_func, produce vectors of sinking speeds at
% depth. Depending on the function definitions, sink speeds may be constant
% with depth, or have a functional relationship specified by a single
% parameter.
Params.wDOM_func = @(wDOM1, nz) wDOM1 .* ones(nz,1);
Params.wDOM1 = 0; % DOM does not sink so wDOM1 should equal 0 and wDOM is constant
Params.wDOM = [];
switch FixedParams.depthDependentPOMsinkSpeed
    case false
        % POM sink speed is constant over depth
        Params.wPOM_func = @(wPOM1, nz) wPOM1 .* ones(nz,1);
        Params.wPOM1 = 1; % this is set far lower than Ward (2012, 2016), but POM sink speed needs to be this low to fit the data, unless remineralisation rates are made far greater...
        Params.wPOM = [];
    case true
        % Allow reduction of POM sinking speed near surface for low values of
        % wPOM1. At large wPOM1 the sink speed is maximal throughout water column.
        % Setting values of wPOM1 on log-scale will improve numerical optimsation
        % of this parameter (the log-scale is smoother across the relevant range),
        % just exponentiate it in wPOM_func.
        Params.wPOM_func = @(wPOM1, wmin, wmax, z) wmax - (wmax - wmin) .* exp(-exp(wPOM1) .* z); % curvature parameter, min and max sink speeds, and depth
        Params.wPOM1 = log(0.15); % curvature of sinking speed-depth curve (log(1/m))
        Params.wPOM = [];
end

Params.wk_func = @(wDOM, wPOM, DOM_i, POM_i) wDOM * DOM_i + wPOM * POM_i; % OM sinking rate (m / day)
Params.wk = [];


%% Parameter bounds

% Not specifically declaring bounds (commenting them out) results in the
% bounds being set equal to the parameter value, i.e, the parameter becomes
% fixed.

% scalars
Bounds.Tref       = [20, 20];
Bounds.A          = [0.01, 0.2];
Bounds.h          = [5, 15];
Bounds.m2         = [0, 1e-3];
% Bounds.aP         = [7.66e-8, 5.75e-6];
Bounds.aP         = [1e-10, 5.75e-6];
% Bounds.aP         = [0, 0.5];
Bounds.theta      = [3, 5];
Bounds.xi         = [1.5, 5];
% Bounds.k_G        = [0.5, 10];
Bounds.k_G        = [0.5, 30];
Bounds.sigG       = [0.25, 2.5];
Bounds.delta_opt  = [10, 10];
Bounds.Lambda     = [-1.5, -0.5];
Bounds.lambda_max = [0.5, 0.9];
Bounds.wDOM1      = [0, 0];

% Bounds.K_m_coef   = [1, 1];
Bounds.K_m_coef   = [0, 0];

% The POM sinking speed has proven to be an awkward parameter... I think
% that including both PON and POC data within the cost function is
% problematic as these combined data create undesirable local minama. I
% think that excluding the POC data from the cost function could be useful,
% or perhaps just fixing wPOM1 to constant value...
switch FixedParams.depthDependentPOMsinkSpeed
    case false
        Bounds.wPOM1 = [0.5, FixedParams.maxSinkSpeed_POM];
    case true
        x = log(2 * (FixedParams.maxSinkSpeed_POM - FixedParams.minSinkSpeed_POM) ./ ...
            FixedParams.maxSinkSpeed_POM);
        Bounds.wPOM1 = log([1 / 30, 1 / 1] .* x); % assume 30m is max depth at which wPOM <= max(wPOM)
end
Bounds.rDOC       = [0.005, 0.06];
Bounds.rDON       = [0.005, 0.06];
Bounds.rPOC       = [0.01, 0.12];
Bounds.rPON       = [0.01, 0.12];

% size dependent

% Quota bounds from Maranon (2013).
% (Bounds derived from the literature are not guarenteed to be consistent 
% with the model constraints, so some inputted values are automatically
% corrected if required).
Bounds.Q_C_a = [Params.Q_C_a, Params.Q_C_a]; % doesn't vary... could maybe go in fixed parameters
Bounds.Q_C_b = [Params.Q_C_b, Params.Q_C_b];

Bounds.Qmin_QC_a = (1/14) * 1e-9 / Params.Q_C_a .* [10^-1.78, 10^-1.26];
Bounds.Qmin_QC_a = max(0, Bounds.Qmin_QC_a); % required Qmin_QC_a > 0

Bounds.Qmin_QC_b = -Params.Q_C_b + [0.77, 0.92];
Bounds.Qmin_QC_b = max(Bounds.Qmin_QC_b, -Params.Q_C_b); % required Qmin_QC_b > - Q_C_b
% Bounds.Qmin_QC_b(2) = 2.5 * Bounds.Qmin_QC_b(2); % Param sensitivity analysis suggest extending upper bound could be useful...

Bounds.Qmax_delQ_a = [10^(-1.78 - (-0.99)), 10^(-1.26 - (-1.35))];
Bounds.Qmax_delQ_a = min(1, max(0, Bounds.Qmax_delQ_a)); % required 0 < Qmax_delQ_a < 1

Bounds.Qmax_delQ_b = [0.77 - 0.96, 0.92 - 0.83];
Bounds.Qmax_delQ_b = min(0, Bounds.Qmax_delQ_b); % required Qmax_delQ_b < 0
% Bounds.Qmax_delQ_b = min(-1e-3, Bounds.Qmax_delQ_b); % required Qmax_delQ_b < 0
% Bounds.Qmax_delQ_b = sort(log(-Bounds.Qmax_delQ_b)); % estimate on log negative scale

Bounds.Vmax_QC_a = 24 / 14 * 1e-9 / Params.Q_C_a * [10^-3.18, 10^-2.78]; % Vmax bounds from Maranon (2013)
Bounds.Vmax_QC_b = -Params.Q_C_b + [0.89, 1.06];
% Bounds.Vmax_QC_b(2) = 3 * Bounds.Vmax_QC_b(2); % param sensitivity suggests extending upper bound could be useful...


% Bounds.aN_QC_a = Bounds.Vmax_QC_a ./ [10^-0.44, 10^-1.2]; % N affinity bounds from Edwards et al. (2015)
% Bounds.aN_QC_b = Bounds.Vmax_QC_b -[0.45, 0.24];

Bounds.aN_QC_a = 1e-3 .* 10 .^ [-9, -7.4] / Params.Q_C_a; % N affinity bounds from Edwards et al. (2015)
Bounds.aN_QC_b = -Params.Q_C_b + [0.58, 0.98];

% % try these more restrictive bounds for affinity -- I think the above provided too much freedom...
% Bounds.aN_QC_a = Params.Vmax_QC_a ./ [10^-0.44, 10^-1.2]; % N affinity bounds from Edwards et al. (2015)
% Bounds.aN_QC_a(2) = 3 * Bounds.aN_QC_a(2); % param sensitivity suggests extending upper bound could be useful...
% Bounds.aN_QC_b = Params.Vmax_QC_b -[0.45, 0.24];
% Bounds.aN_QC_b(2) = -0.025; % param sensitivity suggests extending upper bound could be useful (although size dependency of N affinity should be significant!)
% Bounds.aN_QC_b = sort(log(-Bounds.aN_QC_b)); % estimate on negative log scale


% Bounds.pmax_a = [1.8, 24];  % pmax bounds guessed from mu_inf CIs given in Ward (2017)
% Bounds.pmax_a = [5, 100];
Bounds.pmax_a = [0.5, 5];
Bounds.pmax_b = [-0.5, 0];
% Bounds.pmax_b = [-0.5, -1e-2];
% Bounds.pmax_b = sort(log(-Bounds.pmax_b)); % estimate on negative log scale

% Bounds.Gmax_a = [0, 100];
Bounds.Gmax_a = [5, 35];
Bounds.Gmax_b = [-0.5, 0];
% Bounds.Gmax_b = [-0.5, -1e-2];
% Bounds.Gmax_b = sort(log(-Bounds.Gmax_b)); % esimated on negative log scale

% Bounds.k_G_a = [0, 10];
% Bounds.k_G_b = [0, 1];

Bounds.m_a = [2 .* FixedParams.m_min, 0.1]; % if m_a=m_min then m_b becomes totally irrelevant => set lower bound of m_a > m_min
Bounds.m_b = [-1, 0]; % negativity ensures that mortality rate decreases with size
% Bounds.m_b = [-1, -1e-2]; % negativity ensures that mortality rate decreases with size
% Bounds.m_b = sort(log(-Bounds.m_b)); % estimate on negative log scale

Bounds.beta1 = [0.5, 1];
Bounds.beta2 = [0, 0.9];
Bounds.beta3 = [min(log10(FixedParams.PPsize)), max(log10(FixedParams.PPsize))];

% These bounds ensure that wp never exceeds maxSinkSpeed_P and that
% wp(i)<wp(j) for cell sizes j>i.
Bounds.wp_a = [0, FixedParams.PPsize(end) .^ (-2/3) .* FixedParams.maxSinkSpeed_P];
Bounds.wp_b = [0, log(FixedParams.maxSinkSpeed_P ./ Bounds.wp_a(2)) ./ log(max(FixedParams.PPsize))];


%% Tidy up

[Params, Bounds] = paramCheck(Params, Bounds);


end


%% auxiliarly functions
% Sphere diameter & volume conversions
function vol = d2vol(d)
vol = 1 ./ 6 .* pi .* d .^ 3;
end
% d2vol = @(z) 1 ./ 6 .* pi .* z .^ 3;

% function d = vol2d(vol)
% d = 2 .* (vol .* (3 ./ 4 ./ pi)) .^ (1/3);
% end

function y = powerFunction(a,b,x)
y = a .* x .^ b;
end
% powerFunction = @(a,b,x) a .* x .^ b;
