function [FixedParams, Params] = initialiseParameters(Forc, Data)

% Choose initial values for fixed and variable parameters

%% Fixed Parameters

fields = fieldnames(Forc);
nyears = length(fields);
years = nan(1,nyears);
for i = 1:nyears
    [yy,~] = datevec(Forc.(fields{i}).t);
    years(i) = yy(1);
end

FixedParams.years = years;
FixedParams.nt = length(Forc.(['y' num2str(years(1))]).t); % number of time steps per trajectory - must be fixed over all years

%~~~~~~~~~~~~~~~~~~~~~~~~
% Depth layer arrangement
%~~~~~~~~~~~~~~~~~~~~~~~~
FixedParams.nz = 15;                                                  % number of modelled depth layers
FixedParams.Htot = 150;                                               % total modelled depth
dzmin = 5;                                                            % minimum layer width (dzmin <= Htot/(nz-1))
dzmax = 2 * FixedParams.Htot / FixedParams.nz - dzmin;                % maximum layer width
FixedParams.zw = [0; -cumsum(linspace(dzmin,dzmax,FixedParams.nz))']; % depth of layer edges
FixedParams.zwidth = FixedParams.zw(1:end-1) - FixedParams.zw(2:end); % widths of depth layers
FixedParams.z = 0.5*(FixedParams.zw(1:end-1)+FixedParams.zw(2:end));  % midpoints of depth layers
FixedParams.delz = abs(diff(FixedParams.z));                          % distance between centres of adjacent depth layers

%~~~~~~~~~~~~~~~~
% State variables
%~~~~~~~~~~~~~~~~
% Inorganic nutrients - only nitrate is modelled
FixedParams.IN_nut = {'NO3'};
FixedParams.nIN = length(FixedParams.IN_nut);
% Plankton
FixedParams.PP_nut = {'C','N','Chl'}; % model phytoplankton carbon, nitrogen and chlorophyll
FixedParams.nPP_nut = length(FixedParams.PP_nut);

PPdia = unique(Data.size.size);
PPsize = 4/3*pi*(PPdia ./ 2) .^ 3;
FixedParams.nPP_size = length(PPdia); % number of phytoplankton size classes
FixedParams.PPdia = PPdia; % cell diameter [mu m]
FixedParams.PPsize = PPsize; % cell volume [mu m^3]

% FixedParams.nPP_size = 6; % number of phytoplankton size classes
FixedParams.nPP = FixedParams.nPP_size * FixedParams.nPP_nut;
FixedParams.nZP = 1; % number of zooplankton classes
% Phytoplanton sizes - smallest diameter is 0.5 mu m, volumes of successive size 
% classes increase by factors of 32 (equally spaced on log-scale).
% PPdia = 0.5; % cell diameter [mu m]
% PPsize = 4/3*pi*(PPdia/2)^3; % cell volume [mu m^3]
% PPsize(2:FixedParams.nPP_size) = 32 .^ (1:FixedParams.nPP_size-1) .* PPsize(1);
% PPsize = PPsize(:);
% FixedParams.PPsize = PPsize;
% FixedParams.PPdia = 2 .* (3 .* FixedParams.PPsize ./ (4*pi)) .^ (1/3);
FixedParams.diatoms = FixedParams.PPdia >= 10;                            % assume large phytoplankton are diatoms - only needed to split SINMOD output over classes during state variable initialisation
FixedParams.phytoplankton = [true(1,FixedParams.nPP_size) ... 
    false(1,FixedParams.nZP)]';                                              % index phytoplankton
FixedParams.zooplankton = [false(1,FixedParams.nPP_size) ... 
    true(1,FixedParams.nZP)]';                                               % index zooplankton
% Organic matter
FixedParams.OM_nut = {'C', 'N'};
FixedParams.nOM_nut = length(FixedParams.OM_nut);
FixedParams.OM_type = {'DOM', 'POM'};
FixedParams.nOM_type = length(FixedParams.OM_type);
FixedParams.nOM = FixedParams.nOM_type * FixedParams.nOM_nut;
% All variables
FixedParams.nVar = FixedParams.nIN + FixedParams.nPP + FixedParams.nZP + ...
    FixedParams.nOM;                                                        % number of state variables per depth and trajectory
FixedParams.nEquations = FixedParams.nVar * FixedParams.nz;                 % total number of ODEs

% Indices to extract inorganic nutrients, plankton and organic
% matter from the array containing all variables
FixedParams = createIndexes(FixedParams);

FixedParams.attSW = 0.04; % light attenuation by sea water (1 / m)
FixedParams.attP = 0.0149;  % light attenuation by chlorophyll (m^2 / mg Chl), Krause-Jensen & Sand-Jensen, Limnol. Oceanogr. 43(3):396-407(1998)

FixedParams.POM_is_lost = true; % is POM lost from the system by sinking below bottom modelled depth layer


%% Variable Parameters

% List names of all paramters that may be varied
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
    'wk'
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
    'beta1'
    'beta2'
    'beta3'
    };


%~~~~~~~~~
% Plankton
%~~~~~~~~~

% Size-dependent
% carbon quota [mmol C / cell], Maranon et al. (2013)
Params.Q_C_a = (1/12) * 1e-9 * 10^-0.69;
Params.Q_C_b = 0.88;
Params.Q_C = volumeDependent(Params.Q_C_a, Params.Q_C_b, PPsize);
% min and max cellular nitrogen quota [mmol N / cell], scaled by C quota [mmol N / mmol C], Maranon et al. (2013)
Params.Qmin_QC_a = (1/14) * 1e-9 * 10^-1.47 / Params.Q_C_a;
Params.Qmin_QC_b = 0.84 - Params.Q_C_b;
Params.Qmin_QC = volumeDependent(Params.Qmin_QC_a, Params.Qmin_QC_b, PPsize);
% max quota is parameterised using the ratio
% Qmax/(Qmax-Qmin) = 1/(1-a*V^b), where 0<a<1 and b<0 to ensure Qmax > Qmin
Params.Qmax_delQ_a = 10^(-1.47+1.26);
Params.Qmax_delQ_b = 0.84 - 0.93;
Params.Qmax_delQ = 1 ./ (1 - volumeDependent(Params.Qmax_delQ_a, ...
    Params.Qmax_delQ_b, PPsize));
% maximum nitrogen uptake rate Vmax [mmol N/cell/day] scaled by C quota, Vmax_over_QC [mmol N / mmol C / day], Maranon et al. (2013)
Params.Vmax_QC_a = 24 / 14 * 1e-9 * 10^-3 / Params.Q_C_a;
Params.Vmax_QC_b = 0.97 - Params.Q_C_b;
Params.Vmax_QC = volumeDependent(Params.Vmax_QC_a, Params.Vmax_QC_b, PPsize);

% cellular affinity for nitrogen scaled by QC [m^3 / mmol C / day], derived using half saturation from Litchmann et al. (2007)
Params.aN_QC_a = Params.Vmax_QC_a / 10^-0.77;
Params.aN_QC_b = Params.Vmax_QC_b -0.27;
Params.aN_QC = volumeDependent(Params.aN_QC_a, Params.aN_QC_b, PPsize);

% Params.aN_over_Qmin_a = 24 * 10^(-3 + 0.77 + 1.47);
% Params.aN_over_Qmin_b = 0.97 - 0.27 - 0.84;
% Params.aN_over_Qmin = volumeDependent(Params.aN_over_Qmin_a, ... 
%     Params.aN_over_Qmin_b, PPsize);

% maximum photosynthetic rate [1/day], values guessed based on mu_inf from Ward et al. (2017)
Params.pmax_a = 6;
Params.pmax_b = -0.1;
Params.pmax = volumeDependent(Params.pmax_a, Params.pmax_b, PPsize);
% partitioning of dead matter into DOM and POM
% Params.beta = 0.9 - 0.7 ./ (1 + exp(2.0 - log10(PPsize)));
% Params.beta(FixedParams.nPP_size+1) = Params.beta(FixedParams.nPP_size); % assume beta for zooplankton is equivalent to largest phytoplankton size class

logSize = log10(PPsize);

Params.beta1 = 0.9; % initial values from Ward et al. (2016)
Params.beta2 = 0.2/0.9;
Params.beta3 = 2.0;

expBeta3 = exp(logSize - Params.beta3);
Params.beta = Params.beta1 ./ (1 + expBeta3) .* (1 + Params.beta2 .* expBeta3);

Params.beta(FixedParams.nPP_size+1) = Params.beta(FixedParams.nPP_size); % assume beta for zooplankton is equivalent to largest phytoplankton size class





% Size-independent

% Params.aP = 0.5;          % initial slope of photosynthesis(metabolism)-irradiance curve [mmol N / (mg Chl mu Ein / m^2)], value guessed to produce 'sensible looking' curves given range of irradiances
Params.aP = 7e-3;          % initial slope of photosynthesis(metabolism)-irradiance curve [mmol C / (mg Chl mu Ein / m^2)]
Params.theta = 4.2;         % maximum chl:N ratio [mg Chl / mmol N], Geider et al., 1998
Params.xi = 2.33;         % cost of photosynthesis (mmol C / mmol N)
Params.Tref = 20;         % reference temperature (degrees C)
Params.A = 0.05;          % temperature dependence (unitless)
Params.h = 10;          % curvature on quota uptake limitation
Params.m = 0.05;          % linear plankton mortality (1/day)
Params.k_G = 1;         % half-saturation prey concentration (mmol N / m^3) for grazing uptake
Params.Gmax = 22;          % maximum grazing rate
Params.Lambda = -1;       % prey refuge parameter (unitless)
Params.lambda_max = 0.7;  % maximum prey assimilation efficiency

%~~~~~~~~~~~~~~~
% Organic matter
%~~~~~~~~~~~~~~~
Params.rDON = 0.02;  % dissolved organic nitrogen remineralisation rate (1/day)
Params.rDOC = 0.02;  % dissolved organic carbon remineralisation rate (1/day)
Params.rPON = 0.04;  % particulate organic nitrogen remineralisation rate (1/day)
Params.rPOC = 0.04;  % particulate organic nitrogen remineralisation rate (1/day)
% Params.rPOM = 0.04;  % particulate organic nitrogen remineralisation rate (1/day)
% Params.rDOM = 0.02;  % dissolved organic nitrogen remineralisation rate (1/day)
Params.wk = 10;      % sinking rate of POM (m / day)


%% Functions of parameters

% For efficiency, rescale, reshape and extend dimensions of some parameters
% before using in model. Also calculate any functions of parameters that
% are independent of the state variables
Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;

Params.kN = Params.Vmax_QC ./ Params.aN_QC;

% Params.rDOM = nan(1, FixedParams.nOM_nut);
Params.rOM = nan(FixedParams.nOM_type,1,FixedParams.nOM_nut);
Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_C_index) = Params.rDOC;
Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_N_index) = Params.rDON;
Params.rOM(FixedParams.POM_index,1,FixedParams.OM_C_index) = Params.rPOC;
Params.rOM(FixedParams.POM_index,1,FixedParams.OM_N_index) = Params.rPON;

Params.wk = [0 Params.wk];

% If state variables sink using backwards difference scheme then there's an
% upper limit on integration time steps
dt_max = min(FixedParams.zwidth) ./ Params.wk;
con = true;
tx = 1;
while con
    if any((1 / tx) > dt_max), tx = tx + 1; else, con = false; end
end
FixedParams.dt_max = 1 /tx;


%% Parameter bounds

% Choose bounds to restrict numerical optimisers

% If there is good reason to exclude any particular parameter from the
% optimisation then set its lower/upper bounds equal to its value

% scalars - NEED TO FIND REFERENCES FOR THESE SCALAR BOUNDS...
Params.lowerBound.Tref = 20;        Params.upperBound.Tref = 20;
Params.lowerBound.A = 0.01;         Params.upperBound.A = 0.2;
Params.lowerBound.h = 5;            Params.upperBound.h = 15;
% Params.lowerBound.m = 0.01;         Params.upperBound.m = 0.1;
Params.lowerBound.m = 0.001;         Params.upperBound.m = 0.1;
Params.lowerBound.aP = 0.001;       Params.upperBound.aP = 0.5;
Params.lowerBound.theta = 3;        Params.upperBound.theta = 5;
Params.lowerBound.xi = 3;           Params.upperBound.xi = 5;
Params.lowerBound.Gmax = 1;         Params.upperBound.Gmax = 30;
Params.lowerBound.k_G = 0.01;       Params.upperBound.k_G = 5;
Params.lowerBound.Lambda = -1.5;    Params.upperBound.Lambda = -0.5;
Params.lowerBound.lambda_max = 0.5; Params.upperBound.lambda_max = 0.9;
Params.lowerBound.wk = 10;          Params.upperBound.wk = 10;
% Params.lowerBound.rDOC = 0.01;      Params.upperBound.rDOC = 0.15;
% Params.lowerBound.rDON = 0.01;      Params.upperBound.rDON = 0.15;
% Params.lowerBound.rPOC = 0.01;      Params.upperBound.rPOC = 0.15;
% Params.lowerBound.rPON = 0.01;      Params.upperBound.rPON = 0.15;
Params.lowerBound.rDOC = Params.rDOC * 1/3; Params.upperBound.rDOC = Params.rDOC * 3;
Params.lowerBound.rDON = Params.rDON * 1/3; Params.upperBound.rDON = Params.rDON * 3;
Params.lowerBound.rPOC = Params.rPOC * 1/3; Params.upperBound.rPOC = Params.rPOC * 3;
Params.lowerBound.rPON = Params.rPON * 1/3; Params.upperBound.rPON = Params.rPON * 3;



% size dependent
Params.lowerBound.Qmin_QC_a = (1/14) * 1e-9 * 10^-1.78 / Params.Q_C_a; % Qmin bounds from Maranon (2013)
Params.upperBound.Qmin_QC_a = (1/14) * 1e-9 * 10^-1.26 / Params.Q_C_a;
Params.lowerBound.Qmin_QC_b = 0.77 - Params.Q_C_b;
Params.upperBound.Qmin_QC_b = 0.92 - Params.Q_C_b;

Params.lowerBound.Qmax_delQ_a = 10^(-1.78 - (-0.99)); % Qmax bounds from Maranon (2013)
Params.upperBound.Qmax_delQ_a = 10^(-1.26 - (-1.35));
Params.lowerBound.Qmax_delQ_b = 0.77 - 0.96;
Params.upperBound.Qmax_delQ_b = 0.92 - 0.83;

Params.lowerBound.Vmax_QC_a = 24 / 14 * 1e-9 * 10^-3.18 / Params.Q_C_a; % Vmax bounds from Maranon (2013)
Params.upperBound.Vmax_QC_a = 24 / 14 * 1e-9 * 10^-2.78 / Params.Q_C_a;
Params.lowerBound.Vmax_QC_b = 0.89 - Params.Q_C_b;
Params.upperBound.Vmax_QC_b = 1.06 - Params.Q_C_b;

Params.lowerBound.aN_QC_a = Params.lowerBound.Vmax_QC_a / 10^-0.44; % N affinity bounds from Edwards et al. (2015)
Params.upperBound.aN_QC_a = Params.upperBound.Vmax_QC_a / 10^-1.2;
Params.lowerBound.aN_QC_b = Params.lowerBound.Vmax_QC_b - 0.45;
Params.upperBound.aN_QC_b = Params.upperBound.Vmax_QC_b - 0.24;

% Params.lowerBound.aN_QC_a = Params.lowerBound.Vmax_QC_a / 0.17; % this doesn't, but should, account for CIs of half saturation K intercept a CIs as thsee are not given in Litchman (2013)
% Params.upperBound.aN_QC_a = Params.upperBound.Vmax_QC_a / 0.17;
% Params.lowerBound.aN_QC_b = Params.lowerBound.Vmax_QC_b - 0.36; % N affinity slope bounds from Litchman
% Params.upperBound.aN_QC_b = Params.upperBound.Vmax_QC_b - 0.2;

Params.lowerBound.pmax_a = 1.8; % pmax bounds guessed from mu_inf CIs given in Ward (2017)
Params.upperBound.pmax_a = 24;
Params.lowerBound.pmax_b = -0.7;
Params.upperBound.pmax_b = -0.09;

Params.lowerBound.beta1 = 0.5;
Params.upperBound.beta1 = 1;
Params.lowerBound.beta2 = 0;
Params.upperBound.beta2 = 0.9;
Params.lowerBound.beta3 = min(logSize);
Params.upperBound.beta3 = max(logSize);



end

%%
% power-function of cell volume
function p = volumeDependent(a,b,V)
p = a .* V .^ b;
end



