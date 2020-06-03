function [FixedParams, Params] = initialiseParameters(Forc)

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
FixedParams.INtype = {'NO3'};
FixedParams.nIN = length(FixedParams.INtype);
% Plankton
FixedParams.nPP = 6; % number of phytoplankton size classes
FixedParams.nZP = 1; % number of zooplankton classes
% Phytoplanton sizes - smallest diameter is 0.5 mu m, volumes of successive size 
% classes increase by factors of 32 (equally spaced on log-scale).
PPdia = 0.5;
FixedParams.PPsize = 4/3*pi*(PPdia/2)^3;
FixedParams.PPsize(2:FixedParams.nPP) = 32 .^ (1:FixedParams.nPP-1) .* ...
    FixedParams.PPsize(1);                                                  % cell volumes [mu m^3]
FixedParams.PPsize = FixedParams.PPsize';
FixedParams.PPdia = 2 .* (3 .* FixedParams.PPsize ./ (4*pi)) .^ (1/3);
FixedParams.diatoms = FixedParams.PPsize >= 100;                            % assume large phytoplankton are diatoms - only needed to split SINMOD output over classes during state variable initialisation
FixedParams.phytoplankton = [true(1,FixedParams.nPP) ... 
    false(1,FixedParams.nZP)]';                                              % index phytoplankton
FixedParams.zooplankton = [false(1,FixedParams.nPP) ... 
    true(1,FixedParams.nZP)]';                                               % index zooplankton
% Organic matter - only DOM is explicitly modelled
FixedParams.OMtype = {'DOM'};
FixedParams.nOM = length(FixedParams.OMtype);
% All variables
FixedParams.nVar = FixedParams.nIN + FixedParams.nPP + FixedParams.nZP + ...
    FixedParams.nOM;                                                        % number of state variables per depth and trajectory
FixedParams.nEquations = FixedParams.nVar * FixedParams.nz;                 % total number of ODEs

% Indices to extract inorganic nutrients, plankton and organic
% matter from the array containing all variables
FixedParams = createIndexes(FixedParams);

FixedParams.attSW = 0.04; % light attenuation in sea water
FixedParams.attP = 0.04;  % plankton-specific light attenuation

FixedParams.POM_is_lost = true; % is POM lost from the system by sinking below bottom modelled depth layer

FixedParams.returnExtras = 'auxiliary'; % return extra output from ODEs



%% Variable Parameters

% List names of all paramters that may be varied
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
    'Qmax_over_delQ_a'
    'Qmax_over_delQ_b'
    'Vmax_over_Qmin_a'
    'Vmax_over_Qmin_b'    
    'aN_over_Qmin_a'
    'aN_over_Qmin_b'    
    'pmax_a'
    'pmax_b'
    'beta'
    };


%~~~~~~~~~
% Plankton
%~~~~~~~~~

% Size-dependent

V_PP = FixedParams.PPsize(:); % cell volumes
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
Params.Vmax_over_Qmin_a = 24 * 10^(-3 + 1.47);
Params.Vmax_over_Qmin_b = 0.97 - 0.84;
Params.Vmax_over_Qmin = volumeDependent(Params.Vmax_over_Qmin_a, ... 
    Params.Vmax_over_Qmin_b, V_PP);
% cellular affinity for nitrogen scaled by qmin [m^3 / mmol N / day], derived using half saturation from Litchmann et al. (2007)
Params.aN_over_Qmin_a = 24 * 10^(-3 + 0.77 + 1.26);
Params.aN_over_Qmin_b = 0.97 - 0.27 - 0.84;
Params.aN_over_Qmin = volumeDependent(Params.aN_over_Qmin_a, ... 
    Params.aN_over_Qmin_b, V_PP);
% maximum photosynthetic rate [1/day] at infinite quota, values guessed based on mu_inf from Ward et al. (2017)
Params.pmax_a = 100;
% Params.pmax_a = 35;
Params.pmax_b = -0.26;
Params.pmax = volumeDependent(Params.pmax_a, Params.pmax_b, V_PP);
% partitioning of dead matter into DOM and POM
Params.beta = 0.9 - 0.7 ./ (1 + exp(2.0 - log10(V_PP)));
Params.beta(FixedParams.nPP+1) = Params.beta(FixedParams.nPP); % assume beta for zooplankton is equivalent to largest phytoplankton size class

% Size-independent

Params.aP = 0.5;          % initial slope of photosynthesis(metabolism)-irradiance curve [m^2 / mu Ein], value guessed to produce 'sensible looking' curves given range of irradiances
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

%%
% power-function of cell volume
function p = volumeDependent(a,b,V)
p = a .* V .^ b;
end



