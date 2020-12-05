function [FixedParams, Params, Bounds] = defaultParameters(varargin)

% Set default parameter values.
% Some variables left empty will be assigned values in initialiseParameters.m

if ~isempty(varargin)
    if any(strcmp(varargin, 'Data'))
        % If size spectra data are included in the optional arguments then
        % modelled size classes are chosen to match the data
        Data = varargin{find(strcmp(varargin, 'Data'))+1};
    end
end

%% Fixed parameters

%~~~~~~~~
% Timings
%~~~~~~~~
FixedParams.nt    = []; % number of modelled time steps (determined from forcing data)
FixedParams.years = []; % number of modelled years (determined from forcing and fitting data)

%~~~~~~~~~~~~~
% Depth layers
%~~~~~~~~~~~~~
FixedParams.nz     = 15;                                                   % number of modelled depth layers
FixedParams.Htot   = 150;                                                  % total modelled depth
FixedParams.dzmin  = 5;                                                    % minimum layer width (dzmin <= Htot/(nz-1))
FixedParams.dzmax  = 2 * FixedParams.Htot / FixedParams.nz - ... 
    FixedParams.dzmin;                                                     % maximum layer width
FixedParams.zwidth = linspace(FixedParams.dzmin, FixedParams.dzmax, ... 
    FixedParams.nz)';                                                      % widths of depth layers
FixedParams.zw     = [0; -cumsum(FixedParams.zwidth)];                     % depth of layer edges
FixedParams.zwidth = FixedParams.zw(1:end-1) - FixedParams.zw(2:end);      % widths of depth layers
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
FixedParams.nPP_nut = [];
if ~exist('Data', 'var')
    FixedParams.nPP_size = 9;                                         % number of phytoplankton size classes
    FixedParams.PPdia = 2 .^ (0:FixedParams.nPP_size-1);              % set cell sizes if fitting is not passed as argument
else
    % Automatic if Data is passed as argument
    FixedParams.PPdia = unique(Data.size.size);
    FixedParams.nPP_size = length(FixedParams.PPdia);
end
FixedParams.PPsize        = 4/3 * pi * (FixedParams.PPdia ./ 2) .^ 3; % cell volume [mu m^3]
FixedParams.nPP           = [];                                       % number of phytoplankton classes
FixedParams.nZP           = 1;                                        % number of zooplankton classes
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
FixedParams.attSW        = 0.04;   % light attenuation by sea water (1 / m)
FixedParams.attP         = 0.0149; % light attenuation by chlorophyll (m^2 / mg Chl), Krause-Jensen & Sand-Jensen, Limnol. Oceanogr. 43(3):396-407(1998)
FixedParams.maxSinkSpeed = 10;     % maximum sinking speed (m / day) of POM or cells (this determines max integration time step)
FixedParams.POM_is_lost  = true;   % is POM lost from the system by sinking below bottom modelled depth layer



%% Variable parameters

% List names of all parameters that may be numerically optimised
Params.scalars = {
    'Tref'
    'A'
    'h'
    'm'
    'aP'
    'theta'
    'xi'
    'Gmax'
    'k_G'
    'Lambda'
    'lambda_max'
    'wDOM'
    'wPOM'
    'rDON'
    'rDOC'
    'rPON'
    'rPOC'
    };
Params.sizeDependent = {
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
    'wp_a'
    'wp_b'
    'beta1'
    'beta2'
    'beta3'
    };


%~~~~~~~~~
% Plankton
%~~~~~~~~~

% Size-dependent
% (mostly power functions y=aV^b, formulas given in initialiseParameters.m)

% carbon quota [mmol C / cell], Maranon et al. (2013)
Params.Q_C_a = (1/12) * 1e-9 * 10^-0.69;
Params.Q_C_b = 0.88;
Params.Q_C = [];

% min and max cellular nitrogen quota [mmol N / cell], scaled by C quota [mmol N / mmol C], Maranon et al. (2013)
Params.Qmin_QC_a = (1/14) * 1e-9 * 10^-1.47 / Params.Q_C_a;
Params.Qmin_QC_b = 0.84 - Params.Q_C_b;
Params.Qmin_QC = [];

% max quota
Params.Qmax_delQ_a = 10^(-1.47+1.26);
Params.Qmax_delQ_b = 0.84 - 0.93;
Params.Qmax_delQ = [];

% maximum nitrogen uptake rate Vmax [mmol N/cell/day] scaled by C quota, Vmax_over_QC [mmol N / mmol C / day], Maranon et al. (2013)
Params.Vmax_QC_a = 24 / 14 * 1e-9 * 10^-3 / Params.Q_C_a;
Params.Vmax_QC_b = 0.97 - Params.Q_C_b;
Params.Vmax_QC = [];

% cellular affinity for nitrogen scaled by QC [m^3 / mmol C / day], derived using half saturation from Litchmann et al. (2007)
Params.aN_QC_a = Params.Vmax_QC_a / 10^-0.77;
Params.aN_QC_b = Params.Vmax_QC_b -0.27;
Params.aN_QC = [];

% maximum photosynthetic rate [1/day], values guessed based on mu_inf from Ward et al. (2017)
Params.pmax_a = 6;
Params.pmax_b = -0.1;
Params.pmax = [];

% sinking plankton
Params.wp_a = 5e-3; % intercept = 0 => no sinking at any size
Params.wp_b = 0.1;
Params.wp = [];

% partitioning of dead matter into DOM and POM,  initial values from Ward et al. (2016)
Params.beta1 = 0.9;
Params.beta2 = 0.2/0.9;
Params.beta3 = 2.0;
Params.beta = [];

% Size-independent
Params.aP = 7e-3;           % initial slope of photosynthesis(metabolism)-irradiance curve [mmol C / (mg Chl mu Ein / m^2)]
Params.theta = 4.2;         % maximum chl:N ratio [mg Chl / mmol N], Geider et al., 1998
Params.xi = 2.33;           % cost of photosynthesis (mmol C / mmol N)
Params.Tref = 20;           % reference temperature (degrees C)
Params.A = 0.05;            % temperature dependence (unitless)
Params.h = 10;              % curvature on quota uptake limitation
Params.m = 0.05;            % linear plankton mortality (1/day)
Params.k_G = 1;             % half-saturation prey concentration (mmol N / m^3) for grazing uptake
Params.Gmax = 22;           % maximum grazing rate
Params.Lambda = -1;         % prey refuge parameter (unitless)
Params.lambda_max = 0.7;    % maximum prey assimilation efficiency

%~~~~~~~~~~~~~~~
% Organic matter
%~~~~~~~~~~~~~~~
Params.rDON = 0.02;         % dissolved organic nitrogen remineralisation rate (1/day)
Params.rDOC = 0.02;         % dissolved organic carbon remineralisation rate (1/day)
Params.rPON = 0.04;         % particulate organic nitrogen remineralisation rate (1/day)
Params.rPOC = 0.04;         % particulate organic nitrogen remineralisation rate (1/day)
Params.wDOM = 0;            % POM sinking rate (m / day)
Params.wPOM = 2;


%% Parameter bounds

% Not specifically declaring bounds (commenting them out) results in the
% bounds being set equal to the parameter value, i.e, the parameter becomes
% fixed.

% scalars
Bounds.Tref       = [20, 20];
Bounds.A          = [0.01, 0.2];
Bounds.h          = [5, 15];
Bounds.m          = [0.005, 0.1];
Bounds.aP         = [0.001, 0.5];
Bounds.theta      = [3, 5];
Bounds.xi         = [1.5, 5];
Bounds.Gmax       = [1, 30];
Bounds.k_G        = [0.01, 5];
Bounds.Lambda     = [-1.5, -0.5];
Bounds.lambda_max = [0.5, 0.9];
Bounds.wDOM       = [0, 0];
Bounds.wPOM       = [0, FixedParams.maxSinkSpeed];
Bounds.rDOC       = [0.005, 0.06];
Bounds.rDON       = [0.005, 0.06];
Bounds.rPOC       = [0.005, 0.06];
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

Bounds.Qmax_delQ_a = [10^(-1.78 - (-0.99)), 10^(-1.26 - (-1.35))];
Bounds.Qmax_delQ_a = min(1, max(0, Bounds.Qmax_delQ_a)); % required 0 < Qmax_delQ_a < 1

Bounds.Qmax_delQ_b = [0.77 - 0.96, 0.92 - 0.83];
Bounds.Qmax_delQ_b = min(0, Bounds.Qmax_delQ_b); % required Qmax_delQ_b < 0


Bounds.Vmax_QC_a = 24 / 14 * 1e-9 / Params.Q_C_a * [10^-3.18, 10^-2.78]; % Vmax bounds from Maranon (2013)
Bounds.Vmax_QC_b = -Params.Q_C_b + [0.89, 1.06];

Bounds.aN_QC_a = Bounds.Vmax_QC_a ./ [10^-0.44, 10^-1.2]; % N affinity bounds from Edwards et al. (2015)
Bounds.aN_QC_b = Bounds.Vmax_QC_b -[0.45, 0.24];

Bounds.pmax_a = [1.8, 24];  % pmax bounds guessed from mu_inf CIs given in Ward (2017)
Bounds.pmax_b = [-0.7, -0.09];

Bounds.beta1 = [0.5, 1];
Bounds.beta2 = [0, 0.9];
Bounds.beta3 = [min(log10(FixedParams.PPsize)), max(log10(FixedParams.PPsize))];

Bounds.wp_a = min(FixedParams.zwidth) .* [0.02 ./ 365, 0.05 ./ 1]; % wp_a_range restricts sinking speed of V=1 cells. Lower - some minimal fraction of narrowest depth layer must sink away during 1 year: Upper - some maximum fraction of narrowest depth that sink away in a single day
% r = 2; % choose r>1. assume sinking speed of largest modelled cells is at least r * that of the smallest
% Bounds.wp_b = [log(r) ./ (log(max(FixedParams.PPsize)) - log(min(FixedParams.PPsize))), ...
%     log10(FixedParams.maxSinkSpeed ./ Bounds.wp_a(2)) ./ log10(FixedParams.PPsize(2))];
% Bounds.wp_b(2) = min(2/3, Bounds.wp_b(2)); % speed freely sinking spheres is proportional to D^2, so limit wp_b upper bound

% Assume sinking speed of largest modelled cells is at least 'r' that of
% cells of volume = 1, and bound sinking speed of largest modelled cells to
% maxSinkSpeed.
r = 1; % choose r>=1
Bounds.wp_b = [log10(r) ./ log10(max(FixedParams.PPsize)), ...
    (log10(FixedParams.maxSinkSpeed) - log10(Bounds.wp_a(2))) ./ log10(max(FixedParams.PPsize))];
Bounds.wp_b(2) = min(2/3, Bounds.wp_b(2)); % speed freely sinking spheres is proportional to D^2, so limit wp_b upper bound

%% Tidy up

% Have all parameters listed in scalars and sizeDependent been assigned
% values? If not then display warning
allPars = [Params.scalars; Params.sizeDependent];
for i = 1:length(allPars)
    n = allPars{i};
    try p = Params.(n);
    catch, p = nan;
    end
    if any(isnan(p))
        warning(['Value of ' n ' has not been assigned. See defaultParameters.m'])
        Params.(n) = [];
    end
end

% Flag and correct any out-of-bounds parameters. This just catches any
% incompatabilities between values selected for parameters and their bounds.
% Also, if bounds have not been selected then they are set 
for i = 1:length(allPars)
    n = allPars{i};
    p = Params.(n);
    if ~isempty(p)
        try b = Bounds.(n);
        catch, b = nan;
        end
        if ~any(isnan(b)) % if bounds have been set
            if p < b(1)
                warning(['Default ' n ' value is less than its lower bound, so ' ...
                    n ' has been reset to its lower bound. Choose consistent values in defaultParameters.m'])
                Params.(n) = b(1);
            end
            if p > b(2)
                warning(['Default ' n ' value is greater than its upper bound, so ' ...
                    n ' has been reset to its upper bound. Choose consistent values in defaultParameters.m'])
                Params.(n) = b(2);
            end
        else % if bounds have not been set then parameter is deemed to be fixed
            warning(['Parameter ' n ' declared without bounding values. Bounds now set equal to ' n ' value. See defaultParameters.m'])
            Bounds.(n) = [p, p];
        end
    end
end



