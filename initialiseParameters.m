function [FixedParams, Params] = initialiseParameters(Forc, Data, varargin)

% Load default parameter set and output initialised FixedParams and Params
% structs

%% Load default initial parameters
parFile = [];
if ~isempty(varargin)
    if any(strcmp(varargin, 'load'))
        parFile = varargin{find(strcmp(varargin, 'load'))+1};
    end
end

[FixedParams, Params, Bounds] = defaultParameters('Data', Data, 'load', parFile);

%% Fixed Parameters

%~~~~~~~~~~~~~~~~~~~~~~~~~~
% Timings from forcing data
%~~~~~~~~~~~~~~~~~~~~~~~~~~
fields = fieldnames(Forc);
nf = length(fields);
FixedParams.years = nan(1, nf);
FixedParams.nt = nan(1, nf);
for i = 1:length(fields)
    [y,~] = datevec(Forc.(fields{i}).t);
    FixedParams.years(i) = unique(y);
    FixedParams.nt(i) = length(Forc.(['y' num2str(FixedParams.years(i))]).t);
    if i == length(fields) && length(unique(FixedParams.nt)) > 1
        warning('trajectories should cover same time span each year');
    end
    FixedParams.nt = max(FixedParams.nt);
end

%~~~~~~~~~~~~~
% Depth layers
%~~~~~~~~~~~~~
if isempty(FixedParams.zw)
    FixedParams.zw = [0; -cumsum(linspace(FixedParams.dzmin, ...
        FixedParams.dzmax, FixedParams.nz))'];                               % depth of layer edges
end
if isempty(FixedParams.zwidth)
    FixedParams.zwidth = FixedParams.zw(1:end-1) - FixedParams.zw(2:end);    % widths of depth layers
end
if isempty(FixedParams.z)
    FixedParams.z = 0.5*(FixedParams.zw(1:end-1)+FixedParams.zw(2:end));     % midpoints of depth layers
end
if isempty(FixedParams.delz)
    FixedParams.delz = abs(diff(FixedParams.z));                             % distance between centres of adjacent depth layers
end

%~~~~~~~~~~~~~~~~
% State variables
%~~~~~~~~~~~~~~~~

% Inorganic nutrient
if isempty(FixedParams.nIN)
    FixedParams.nIN = length(FixedParams.IN_nut);  % number of nutrient types
end

% Plankton
if isempty(FixedParams.nPP_nut)
    FixedParams.nPP_nut = length(FixedParams.PP_nut);              % number of nutrient types
end
if isempty(FixedParams.PPdia)
    FixedParams.PPdia = 2 .^ (0:FixedParams.nPP_size-1);           % cell diameter
end
if isempty(FixedParams.PPsize)
    FixedParams.PPsize = 4/3 * pi * (FixedParams.PPdia ./ 2) .^ 3; % cell volume [mu m^3]
end
if isempty(FixedParams.nPP)
    FixedParams.nPP = FixedParams.nPP_size * FixedParams.nPP_nut;  % number of phytoplankton classes
end
if isempty(FixedParams.diatoms)
    FixedParams.diatoms = FixedParams.PPdia >= 10;                 % assume large phytoplankton are diatoms - only needed to split SINMOD output over classes during state variable initialisation
end
if isempty(FixedParams.phytoplankton)
    FixedParams.phytoplankton = [true(1,FixedParams.nPP_size) ... 
        false(1,FixedParams.nZP)]';                                % index phytoplankton
end
if isempty(FixedParams.zooplankton)
    FixedParams.zooplankton = [false(1,FixedParams.nPP_size) ... 
    true(1,FixedParams.nZP)]';                                     % index zooplankton 
end

% Organic matter
if isempty(FixedParams.nOM_nut)
    FixedParams.nOM_nut = length(FixedParams.OM_nut);             % number of nutrient types
end
if isempty(FixedParams.nOM_type)
    FixedParams.nOM_type = length(FixedParams.OM_type);           % number OM types (DOM and POM)
end
if isempty(FixedParams.nOM)
    FixedParams.nOM = FixedParams.nOM_type * FixedParams.nOM_nut; % total number of OM variables
end

% All variables
if isempty(FixedParams.nVar)
    FixedParams.nVar = FixedParams.nIN + FixedParams.nPP + ... 
        FixedParams.nZP + FixedParams.nOM;                        % number of state variables per depth and trajectory
end
if isempty(FixedParams.nEquations)
    FixedParams.nEquations = FixedParams.nVar * FixedParams.nz;   % total number of ODEs
end


% Indices to extract inorganic nutrients, plankton and organic
% matter from the array containing all variables
FixedParams = createIndexes(FixedParams);

% If state variables sink using backwards difference scheme then there's an
% upper limit on integration time steps
dt_max = min(FixedParams.zwidth) ./ ([FixedParams.maxSinkSpeed_POM, ... 
    FixedParams.maxSinkSpeed_P]);
con = true;
tx = 1;
while con
    % set max integration time step as a fraction of a day
    if any((1 / tx) > dt_max), tx = tx + 1; else, con = false; end
end
FixedParams.dt_max = 1 /tx;


%% Variable Parameters

% Size-dependent
Vol = FixedParams.PPsize;
if isempty(Params.Q_C)
    Params.Q_C = powerFunction(Params.Q_C_a, Params.Q_C_b, Vol);
end
if isempty(Params.Qmin_QC)
    Params.Qmin_QC = powerFunction(Params.Qmin_QC_a, Params.Qmin_QC_b, Vol);
end
if isempty(Params.Qmax_delQ)
    Params.Qmax_delQ = 1 ./ (1 - powerFunction(Params.Qmax_delQ_a, ... 
        Params.Qmax_delQ_b, Vol));
end
if isempty(Params.Vmax_QC)
    Params.Vmax_QC = powerFunction(Params.Vmax_QC_a, Params.Vmax_QC_b, Vol);
end
if isempty(Params.aN_QC)
    Params.aN_QC = powerFunction(Params.aN_QC_a, Params.aN_QC_b, Vol);
end
if isempty(Params.pmax)
    Params.pmax = powerFunction(Params.pmax_a, Params.pmax_b, Vol);
end

if isempty(Params.Gmax)
    Params.Gmax = powerFunction(Params.Gmax_a, Params.Gmax_b, Vol);
end

if isempty(Params.wp)
    Params.wp = powerFunction(Params.wp_a, Params.wp_b, Vol, ... 
        'UpperBound', FixedParams.maxSinkSpeed_P, 'Transpose', true);
end
if isempty(Params.beta)
    Params.beta = doubleLogisticFunction(Params.beta1, Params.beta2, ...
        Params.beta3, log10(Vol)); % flexible 3-parameter double logistic function
end


% Functions of parameters
Params.Qmax_QC = Params.Qmin_QC ./ (1 - 1 ./ Params.Qmax_delQ);
Params.delQ_QC = Params.Qmax_QC - Params.Qmin_QC;

Params.kN = Params.Vmax_QC ./ Params.aN_QC;

Params.beta(FixedParams.nPP_size+1) = Params.beta(FixedParams.nPP_size); % assume beta for zooplankton is equivalent to largest phytoplankton size class

Params.rOM = nan(FixedParams.nOM_type,1,FixedParams.nOM_nut);
Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_C_index) = Params.rDOC;
Params.rOM(FixedParams.DOM_index,1,FixedParams.OM_N_index) = Params.rDON;
Params.rOM(FixedParams.POM_index,1,FixedParams.OM_C_index) = Params.rPOC;
Params.rOM(FixedParams.POM_index,1,FixedParams.OM_N_index) = Params.rPON;

Params.wk = nan(1,2);
Params.wk(FixedParams.DOM_index) = Params.wDOM;
Params.wk(FixedParams.POM_index) = Params.wPOM;




%% Parameter bounds
if exist('Bounds', 'var')
    Params.bounds = Bounds;
end




end


function y = powerFunction(a, b, x, varargin)
y = a .* x .^ b;
if ~isempty(varargin)
    tr = strcmp(varargin, 'Transpose');
    if any(tr) && varargin{find(tr)+1}, y = y'; end
    ub = strcmp(varargin, 'UpperBound');
    if any(ub)
        cap = varargin{find(ub)+1};
        y(y>cap) = cap;
    end
end
end

function y = doubleLogisticFunction(a, b, c, x)
u = exp(x - c);
y = a ./ (1 + u) .* (1 + b .* u);
end






