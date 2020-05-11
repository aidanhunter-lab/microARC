function [Params, FixedParams] = initialiseParameters(FixedParams)


%% Parameter names
Params.scalars = {
    'Tref'
    'A'
    'm'
    'aP'
    'Gmax'
    'k_G'
    'Lambda'
    'lambda_max'
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
    'beta'
    };


%% Parameter values

%~~~~~~~~~
% Plankton
%~~~~~~~~~

% Size-dependent

V_PP = FixedParams.PPsize; % cell volumes

% minimum and maximum cellular nitrogen quota [mmol N / cell], values from Maranon et al. (2013)
Params.Qmin_a = (1/14) * 1e-9 * 10^-1.47;
Params.Qmin_b = 0.84;
Params.Qmin = volumeDependent(Params.Qmin_a, Params.Qmin_b, V_PP);

% maximumm quota is parameterised using the ratio
% Qmax/(Qmax-Qmin) = 1/(1-a*V^b), where 0<a<1 and b<0
Params.Qmax_over_delQ_a = 10^(-1.47+1.26);
Params.Qmax_over_delQ_b = 0.84 - 0.93;
Params.Qmax_over_delQ = 1 ./ (1 - volumeDependent(Params.Qmax_over_delQ_a, ...
    Params.Qmax_over_delQ_b, V_PP));

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
Params.pmax = volumeDependent(Params.pmax_a, Params.pmax_b, V_PP);

% partitioning of dead matter into DOM and POM
Params.beta = 0.9 - 0.7 ./ (1 + exp(2.0 - log10(V_PP)));
Params.beta(FixedParams.nPP+1) = Params.beta(FixedParams.nPP); % assume beta for zooplankton is equivalent to largest phytoplankton size class


% Size-independent

% initial slope of photosynthesis(metabolism)-irradiance curve [m^2 / mu Ein], value guessed to produce 'sensible looking' curves given range of irradiances
Params.aP = 0.5;

Params.Tref = 20;         % reference temperature (degrees C)
Params.A = 0.05;          % temperature dependence (unitless)
Params.m = 0.05;          % linear plankton mortality (1/day)
Params.k_G = 0.1;         % half-saturation prey concentration (mmol N / m^3) for grazing uptake
Params.Gmax = 1;          % maximum grazing rate
Params.Lambda = -1;       % prey refuge parameter (unitless)
Params.lambda_max = 0.7;  % maximum prey assimilation efficiency


%~~~~~~~~~~~~~~~
% Organic matter
%~~~~~~~~~~~~~~~

Params.rPOM = 0.04;  % particulate organic nitrogen remineralisation rate (1/day)
Params.rDOM = 0.02;  % dissolved organic nitrogen remineralisation rate (1/day)
Params.wk = 10;      % sinking rate of POM (m / day)


%% Functions of parameters

% For efficiency, reshape and extend dimensions of some parameters before
% using in model. Also calculate functions whose arguments only involve
% parameters and not state variables.


Params.Qmin = Params.Qmin';
Params.Qmax_over_delQ = Params.Qmax_over_delQ';
Params.Vmax_over_qmin = Params.Vmax_over_qmin';
Params.aN_over_qmin = Params.aN_over_qmin';
Params.pmax = Params.pmax';
Params.beta = Params.beta';

% POM sinking and remineralisation matrix
sinkTime = FixedParams.delz ./ Params.wk;         % time particles take to sink from center of one depth layer to the next
r_x_sinkTime = Params.rPOM .* sinkTime;
nz = FixedParams.nz;
dzm = FixedParams.zwidth(2:nz) ./ ...
    (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));
dzp = FixedParams.zwidth(1:nz-1) ./ ...
    (FixedParams.zwidth(1:nz-1) + FixedParams.zwidth(2:nz));

POM_to_IN_array = zeros(nz, nz); % (i,j) lower-tri matrix of proportions of POM remineralised while sinking from layer j to i
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

Params.POM_to_IN = sparse(POM_to_IN); % using sparse matrix is increasingly useful the more nutrient s are modelled...

end


function p = volumeDependent(a,b,V)
p = a .* V .^ b;
end