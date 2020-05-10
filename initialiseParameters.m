function [Params, FixedParams] = initialiseParameters_MonodModelSinglePredator_properParams(FixedParams)

% Store names of parameters whose values may need updated

Params.scalars = {
    'Tref'
    'A'
%     'alpha'
%     'xi'
    'm'
    'aP'
    'Gmax'
    'k_G'
    'Lambda'
%     'sigma'
%     'delta_opt'
    'lambda_max'
%     'h'
    'wk'
    'rPOM'
    'rDOM'
    };

Params.sizeDependent = {
    'Qmin_a'
    'Qmin_b'
    'Qmax_over_delq_a'
    'Qmax_over_delq_b'
    'Vmax_over_qmin_a'
    'Vmax_over_qmin_b'    
    'aN_over_qmin_a'
    'aN_over_qmin_b'    
    'pmax_a'
    'pmax_b'
%     'Gmax_a'
%     'Gmax_b'
    'beta'
    };



% %~~~~~~~~~~~~~~~~~~~~
% % Inorganic nutrients
% %~~~~~~~~~~~~~~~~~~~~
% 
% Params.zeta_NH4 = 2;    % ammonium to nitrite oxidation rate (1/day)
% Params.zeta_NO2 = 0.1;  % nitrite to nitrate oxidation rate (1/day)
% Params.IOx      = 10;   % upper PAR threshold for nitrification (micro Ein / s / m2)

%~~~~~~~~~
% Plankton
%~~~~~~~~~

% Size-dependent

V_PP = FixedParams.PPsize;
% V_ZP = FixedParams.ZPsize;
% V_all = FixedParams.size;

% minimum and maximum cellular nitrogen quota [mmol N / cell], values from Maranon et al. (2013)
Params.Qmin_a = (1/14) * 1e-9 * 10^-1.47;
Params.Qmin_b = 0.84;
Params.Qmin = volumeDependent(Params.Qmin_a, Params.Qmin_b, V_PP); % minimum cell quota (mmol N / cell)

% maximumm quota is parameterised using the ratio
% Qmax/(Qmax-Qmin) = 1/(1-a*V^b), 0<a<1 and b<0
Params.Qmax_over_delQ_a = 10^(-1.47+1.26);
Params.Qmax_over_delQ_b = 0.84 - 0.93;
Params.Qmax_over_delQ = 1 ./ (1 - volumeDependent(Params.Qmax_over_delQ_a, ...
    Params.Qmax_over_delQ_b, V_PP));

% Params.Qmax_a = (1/14) * 1e-9 * 10^-1.26;
% Params.Qmax_b = 0.93;
% Params.Qmax = volumeDependent(Params.Qmax_a, Params.Qmax_b, V_PP); % maximum cell quota (mmol N / cell)


% nitrogen specific maximum uptake rate [1/day], values from Maranon et al. (2013) scaled by qmin
Params.Vmax_over_qmin_a = 24 * 10^(-3 + 1.47);
Params.Vmax_over_qmin_b = 0.97 - 0.84;
Params.Vmax_over_qmin = volumeDependent(Params.Vmax_over_qmin_a, ... 
    Params.Vmax_over_qmin_b, V_PP);

% cellular affinity for nitrogen scaled by qmin [m^3 / mmol N / day], derived using half saturation from Litchmann et al. (2007)
Params.aN_over_qmin_a = 24 * 10^(-3 + 0.77 + 1.26);
Params.aN_over_qmin_b = 0.97 - 0.27 - 0.84;
Params.aN_over_qmin = volumeDependent(Params.aN_over_qmin_a, ... 
    Params.aN_over_qmin_b, V_PP);

% maximum photosynthetic rate [1/day] at infinite quota, values guessed based on mu_inf from Ward et al. (2017)
Params.pmax_a = 100;
% Params.pmax_a = 35;
Params.pmax_b = -0.26;
Params.pmax = volumeDependent(Params.pmax_a, Params.pmax_b, V_PP); % maximum photosynthetic rate (mmol N / cell / day)

% initial slope of photosynthesis(metabolism)-irradiance curve [m^2 / mu Ein], value guessed to produce 'sensible looking' curves given range of irradiances
Params.aP = 0.5;
% Params.aP = 0.05;



Params.beta = 0.9 - 0.7 ./ (1 + exp(2.0 - log10(V_PP))); % partitioning of dead matter into DOM and POM
Params.beta(FixedParams.nPP+1) = Params.beta(FixedParams.nPP); % assume beta for zooplankton is equivalent to largest phytoplankton size class


% Size-independent
Params.Tref       = 20;       % reference temperature (degrees C)
Params.A          = 0.05;     % temperature dependence (unitless)
Params.m         = 0.05;      % linear plankton mortality (1/day)

% Params.alphap = 7.6*10^-6;    % initial gradient in photosynthesis-irradiance curve (mmol N m^2 / cell / muEin)
Params.k_G        = 0.1;        % half-saturation prey concentration (mmol N / m^3) for grazing uptake
Params.Gmax = 1;              % maximum grazing rate
Params.Lambda     = -1;       % prey refuge parameter (unitless)
% Params.sigma      = 2;        % log-scale standard deviation of prey-size preference (unitless)
% Params.delta_opt  = 10;       % optimum predator:prey radius ratio (unitless)
Params.lambda_max = 0.7;      % maximum prey assimilation efficiency
% Params.h          = 0.1;      % shape parameter for nitrogen inorganic nutirnet uptakes


%~~~~~~~~~~~~~~~
% Organic matter
%~~~~~~~~~~~~~~~

% carbon   = strcmp(FixedParams.OMnut, 'C');
% nitrogen = strcmp(FixedParams.OMnut, 'N');
% DOM      = strcmp(FixedParams.OMtype, 'DOM');
% POM      = strcmp(FixedParams.OMtype, 'POM');

% Params.r_POM_C = 0.04;  % particulate organic carbon remineralisation rate (1/day)
Params.rPOM = 0.04;  % particulate organic nitrogen remineralisation rate (1/day)
% Params.r_DOM_C = 0.02;  % dissolved organic carbon remineralisation rate (1/day)
Params.rDOM = 0.02;  % dissolved organic nitrogen remineralisation rate (1/day)

% Params.rDOM(carbon)   = Params.r_DOM_C;
% Params.rDOM(nitrogen) = Params.r_DOM_N;
% Params.rPOM(carbon)   = Params.r_POM_C;
% Params.rPOM(nitrogen) = Params.r_POM_N;

% Params.r(carbon, 1)   = Params.r_DOM_C;
% Params.r(nitrogen, 1) = Params.r_DOM_N;
% Params.r(carbon, 2)   = Params.r_POM_C;
% Params.r(nitrogen, 2) = Params.r_POM_N;

Params.wk = 10;          % sinking rate of POM (m / day)
% Params.wk_POM = 10;    % sinking rate of POM (m / day)
% Params.wk_DOM = 0;     % sinking rate of DOM (m / day)
% Params.wk(POM,1) = Params.wk_POM;
% Params.wk(DOM,1) = Params.wk_DOM;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For efficiency, reshape and extend dimensions of some parameters before
% using in model. Also calculate functions whose arguments only involve
% parameters and not state variables.

% nZP = FixedParams.nZP;

% Params.Qmax_m_Qmin = Params.Qmax - Params.Qmin;                              % nutrient storage capacity is difference between Qmax and Qmin

Params.Qmin = Params.Qmin'; % transpose
Params.Qmax_over_delQ = Params.Qmax_over_delQ';
Params.Vmax_over_qmin = Params.Vmax_over_qmin';
Params.aN_over_qmin = Params.aN_over_qmin';
Params.pmax = Params.pmax';

% Params.Qmin = [Params.Qmin zeros(1, 1)]'; % transpose, and pad with zeros for zooplankton
% Params.Qmax_over_delQ = [Params.Qmax_over_delQ zeros(1, 1)]';
% Params.Vmax_over_qmin = [Params.Vmax_over_qmin zeros(1, 1)]'; % transpose, and pad with zeros for zooplankton
% Params.aN_over_qmin = [Params.aN_over_qmin zeros(1, 1)]';
% Params.pmax = [Params.pmax zeros(1, 1)]';    % transpose, and pad with zeros for zooplankton


Params.beta = Params.beta';


% POM sinking and remineralisation matrix
sinkTime = FixedParams.delz ./ Params.wk;         % time particles take to sink from center of one depth layer to the next
r_x_sinkTime = Params.rPOM .* sinkTime;
nz = FixedParams.nz;
dzm = FixedParams.zwidth(2:nz) ./ ...
    (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
dzp = FixedParams.zwidth(1:nz-1) ./ ...
    (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
% p = zeros(nz-1, nz-1); % proportions of POM remineralised during sinking, at grid points

POM_to_IN_array = zeros(nz, nz); % same as p but at grid centers
% POM_to_IN = zeros(nz * FixedParams.nOM_nut, nz * FixedParams.nOM_nut); % same as p but at grid centers - a block-diagonal matrix
POM_to_IN = []; % create a block-diagonal matrix - only needed when modelling multiple nutrients
for i_nut = 1:FixedParams.nOM
    p = r_x_sinkTime .* tril(cumprod(tril(1 - repmat(r_x_sinkTime, [1 nz-1]), -1) + triu(ones(nz-1,nz-1))));
    POM_to_IN_array(1:nz-1,1:nz-1) = dzp .* p;
    POM_to_IN_array(2:nz,1:nz-1) = POM_to_IN_array(2:nz,1:nz-1) + dzm .* p;
    if ~FixedParams.POM_is_lost
        % If POM does not sink below bottom depth layer then all that
        % remains after sinking is remineralised on bottom layer
        POM_to_IN_array(nz,:) = POM_to_IN_array(nz,:) + (1 - sum(POM_to_IN_array));
    end
    
    POM_to_IN = blkdiag(POM_to_IN, POM_to_IN_array);
end

Params.POM_to_IN = sparse(POM_to_IN);

end


function p = volumeDependent(a,b,V)
p = a .* V .^ b;
end