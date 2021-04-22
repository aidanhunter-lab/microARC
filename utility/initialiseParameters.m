function [FixedParams, Params] = initialiseParameters(Forc, bioModel, varargin)
% Load default parameter set. Output FixedParams and Params structs.

extractVarargin(varargin)

if ~exist('parFile', 'var')
    parFile = [];
end

if ~isempty(parFile)
    useSavedPars = true;
    savedpars = matfile(parFile, 'Writable', true);
    savedpars = savedpars.pars;
else
    useSavedPars = false;
end

if ~exist('nsizes', 'var')
    nsizes = []; % number of plankton size classes
else
    nsizes = int16(eval('nsizes')); % convert to integer type
end

if ~exist('ESDmin', 'var')
    ESDmin = []; % min cell size
end
if ~exist('ESDmax', 'var')
    ESDmax = []; % max cell size
end

% Create FixedParams and Params structs by calling default values
[FixedParams, Params, Bounds] = defaultParameters(bioModel);
% [FixedParams, Params, Bounds] = defaultParameters(bioModel, 'Data', Data, ... 
%     'load', parFile);

% Replace default Params if using saved values
if useSavedPars
    fields = fieldnames(Params);
    for i = 1:height(savedpars)
        if ~any(strcmp(savedpars.Param{i}, fields))
            error(['Loaded parameters mismatch with the default set: ' savedpars.Param{i} ' is not present in the default parameters'])
        else
            Params.(savedpars.Param{i}) = savedpars.Value(i);
        end
    end
end


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
if isempty(FixedParams.nZP_nut)
    FixedParams.nZP_nut = length(FixedParams.ZP_nut);              % number of nutrient types
end


% Set plankton cell sizes using the size-spectra data, if these data were
% passed as an argument in varargin.

if ~isempty(nsizes)
    FixedParams.nPP_size = nsizes;
    if ~isempty(ESDmin) && ~isempty(ESDmax)
        lim = log10([ESDmin; ESDmax]);
        edges = linspace(lim(1), lim(2), nsizes + 1)';
        FixedParams.PPdia_intervals = 10 .^ edges;
        FixedParams.PPdia = 10 .^ (0.5 .* (edges(1:end-1) + edges(2:end)));
    else
        FixedParams.PPdia = 2 .^ (0:nsizes-1);
        FixedParams.PPdia = FixedParams.PPdia(:);
    end
    FixedParams.PPsize = d2vol(FixedParams.PPdia);
end

% Set predator size classes equal to phytoplankton sizes
FixedParams.nZP_size = FixedParams.nPP_size;
FixedParams.ZPdia_intervals = FixedParams.PPdia_intervals;
FixedParams.ZPdia = FixedParams.PPdia;
FixedParams.ZPsize = FixedParams.PPsize;

delta_pred = reshape(FixedParams.ZPsize, [FixedParams.nZP_size, 1]);
delta_prey = [reshape(FixedParams.PPsize, [1 FixedParams.nPP_size]), ... 
    reshape(FixedParams.ZPsize, [1 FixedParams.nZP_size])];
FixedParams.delta = delta_pred ./ delta_prey; % predator:prey size ratios


% if isempty(FixedParams.PPdia)
%     FixedParams.PPdia = 2 .^ (0:FixedParams.nPP_size-1);           % cell diameter
% end
% if isempty(FixedParams.PPsize)
%     FixedParams.PPsize = 4/3 * pi * (FixedParams.PPdia ./ 2) .^ 3; % cell volume [mu m^3]
% end
if isempty(FixedParams.nPP)
    FixedParams.nPP = FixedParams.nPP_size * FixedParams.nPP_nut;  % number of phytoplankton classes
end

% if isempty(FixedParams.ZPdia)
%     FixedParams.ZPdia = 2 .^ (0:FixedParams.nZP_size-1);           % cell diameter
% end
% if isempty(FixedParams.ZPsize)
%     FixedParams.ZPsize = 4/3 * pi * (FixedParams.ZPdia ./ 2) .^ 3; % cell volume [mu m^3]
% end
if isempty(FixedParams.nZP)
    FixedParams.nZP = FixedParams.nZP_size * FixedParams.nZP_nut;  % number of zooplankton classes
end

FixedParams.diaAll = [FixedParams.PPdia; FixedParams.ZPdia];
FixedParams.sizeAll = [FixedParams.PPsize; FixedParams.ZPsize];


if isempty(FixedParams.diatoms)
    FixedParams.diatoms = FixedParams.diaAll >= 10;                 % assume large phytoplankton are diatoms - only needed to split SINMOD output over classes during state variable initialisation
end
if isempty(FixedParams.phytoplankton)
    FixedParams.phytoplankton = [true(1,FixedParams.nPP_size) ... 
        false(1,FixedParams.nZP_size)]';                                % index phytoplankton
end
if isempty(FixedParams.zooplankton)
    FixedParams.zooplankton = [false(1,FixedParams.nPP_size) ... 
    true(1,FixedParams.nZP_size)]';                                     % index zooplankton 
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
Vol = FixedParams.sizeAll;
% VolP = FixedParams.PPsize;
VolZ = FixedParams.ZPsize;

if (isfield(Params, 'Q_C_a') && isfield(Params, 'Q_C_b')) && ...
        ~(isempty(Params.Q_C_a) || isempty(Params.Q_C_b))
    Params.Q_C = Params.Q_C_func(Params.Q_C_a, Params.Q_C_b, Vol);
end

if (isfield(Params, 'Qmin_QC_a') && isfield(Params, 'Qmin_QC_b')) && ...
        ~(isempty(Params.Qmin_QC_a) || isempty(Params.Qmin_QC_b))
    Params.Qmin_QC = Params.Qmin_QC_func(Params.Qmin_QC_a, Params.Qmin_QC_b, Vol);
end

if (isfield(Params, 'Qmax_delQ_a') && isfield(Params, 'Qmax_delQ_b')) && ...
        ~(isempty(Params.Qmax_delQ_a) || isempty(Params.Qmax_delQ_b))
    Params.Qmax_delQ = Params.Qmax_delQ_func(Params.Qmax_delQ_a, Params.Qmax_delQ_b, Vol);
end

if (isfield(Params, 'Vmax_QC_a') && isfield(Params, 'Vmax_QC_b')) && ...
        ~(isempty(Params.Vmax_QC_a) || isempty(Params.Vmax_QC_b))
    Params.Vmax_QC = Params.Vmax_QC_func(Params.Vmax_QC_a, Params.Vmax_QC_b, Vol);
end

if (isfield(Params, 'aN_QC_a') && isfield(Params, 'aN_QC_b')) && ...
        ~(isempty(Params.aN_QC_a) || isempty(Params.aN_QC_b))
    Params.aN_QC = Params.aN_QC_func(Params.aN_QC_a, Params.aN_QC_b, Vol);
end

if (isfield(Params, 'pmax_a') && isfield(Params, 'pmax_b')) && ...
        ~(isempty(Params.pmax_a) || isempty(Params.pmax_b))
    Params.pmax = Params.pmax_func(Params.pmax_a, Params.pmax_b, Vol);
end

if (isfield(Params, 'Gmax_a') && isfield(Params, 'Gmax_b')) && ...
        ~(isempty(Params.Gmax_a) || isempty(Params.Gmax_b))
    switch bioModel
        case 'singlePredatorClass'
            Params.Gmax = Params.Gmax_func(Params.Gmax_a, Params.Gmax_b, Vol)';
        case {'multiplePredatorClasses','mixotrophy'}
            Params.Gmax = Params.Gmax_func(Params.Gmax_a, Params.Gmax_b, VolZ);
    end
end

% if isempty(Params.k_G)
%     Params.k_G = powerFunction(Params.k_G_a, Params.k_G_b, VolZ);
% end

if (isfield(Params, 'm_a') && isfield(Params, 'm_b')) && ...
        ~(isempty(Params.m_a) || isempty(Params.m_b))
    Params.m = Params.m_func(Params.m_a, Params.m_b, FixedParams.m_min, Vol);
end

if (isfield(Params, 'wp_a') && isfield(Params, 'wp_b')) && ...
        ~(isempty(Params.wp_a) || isempty(Params.wp_b))
    Params.wp = Params.wp_func(Params.wp_a, Params.wp_b, Vol);
end

if (isfield(Params, 'beta1') && isfield(Params, 'beta2') && isfield(Params, 'beta3')) && ...
        ~(isempty(Params.beta1) || isempty(Params.beta2) || isempty(Params.beta3))
    Params.beta = Params.beta_func(Params.beta1, Params.beta2, Params.beta3, Vol);
end


% Aggregate some parameters into matrices/arrays

if (isfield(Params, 'wDOM') && isfield(Params, 'wPOM')) && ...
        ~(isempty(Params.wDOM) || isempty(Params.wPOM))
    Params.wk = Params.wk_func(Params.wDOM, Params.wPOM, ... 
        FixedParams.DOM_index, FixedParams.POM_index);
end

if (isfield(Params, 'rDOC') && isfield(Params, 'rDON') && isfield(Params, 'rPOC') && isfield(Params, 'rPON')   ) && ...
        ~(isempty(Params.rDOC) || isempty(Params.rDON) || isempty(Params.rPOC) || isempty(Params.rPON) )
    Params.rOM = Params.rOM_func(Params.rDOC, Params.rDON, Params.rPOC, Params.rPON, ... 
        FixedParams.DOM_index, FixedParams.POM_index, FixedParams.OM_C_index, FixedParams.OM_N_index);
end


% Functions of parameters (not estimated directly)

if (isfield(Params, 'Qmin_QC') && isfield(Params, 'Qmax_delQ')) && ...
        ~(isempty(Params.Qmin_QC) || isempty(Params.Qmax_delQ))
    Params.Qmax_QC = Params.Qmax_QC_func(Params.Qmin_QC, Params.Qmax_delQ);
end

if (isfield(Params, 'Qmin_QC') && isfield(Params, 'Qmax_delQ')) && ...
        ~(isempty(Params.Qmin_QC) || isempty(Params.Qmax_delQ))
    Params.delQ_QC = Params.delQ_QC_func(Params.Qmin_QC, Params.Qmax_QC);
end

if (isfield(Params, 'Vmax_QC') && isfield(Params, 'aN_QC')) && ...
        ~(isempty(Params.Vmax_QC) || isempty(Params.aN_QC))
    Params.kN = Params.kN_func(Params.Vmax_QC, Params.aN_QC);
end

% Params.beta(FixedParams.nPP_size+1) = Params.beta(FixedParams.nPP_size); % assume beta for zooplankton is equivalent to largest phytoplankton size class


%% Tidy up

% Some parameter bounds are dependent upon choice of cell sizes.
% These bounds need updated if chosen cell sizes differ from the defaults.
Bounds.beta3 = [min(log10(FixedParams.PPsize)), max(log10(FixedParams.PPsize))];
Bounds.wp_a = [0, FixedParams.PPsize(end) .^ (-2/3) .* FixedParams.maxSinkSpeed_P];
Bounds.wp_b = [0, log(FixedParams.maxSinkSpeed_P ./ Bounds.wp_a(2)) ./ log(max(FixedParams.PPsize))];

[Params, Bounds] = paramCheck(Params, Bounds);
Params.bounds = Bounds;

end


%% auxiliary functions
% function y = powerFunction(a, b, x, varargin)
% y = a .* x .^ b;
% if ~isempty(varargin)
%     tr = strcmp(varargin, 'Transpose');
%     if any(tr) && varargin{find(tr)+1}, y = y'; end
%     ub = strcmp(varargin, 'UpperBound');
%     if any(ub)
%         cap = varargin{find(ub)+1};
%         y(y>cap) = cap;
%     end
% end
% end
% 
% function y = doubleLogisticFunction(a, b, c, x)
% u = exp(x - c);
% y = a ./ (1 + u) .* (1 + b .* u);
% end

function vol = d2vol(d)
vol = 4 ./ 3 .* pi .* (0.5 .* d) .^ 3;
end

% function d = vol2d(vol)
% d = 2 .* (vol .* (3 ./ 4 ./ pi)) .^ (1/3);
% end

