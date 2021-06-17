function [ESDP, ESDZ, VolP, VolZ, ESDPedges, ESDZedges] = ... 
    createSizeClasses(n, ESDmin, ESDmax, ESDP1, deltaOpt)
% Returns size class grids: cell ESDs and volumes for predators and prey
% set so that size classes are equally spaced on log-scale; that all prey
% classes are optimally grazed by some predator class; and all size classes
% are within set limits.

% n        = desired number of prey size classes -- may be only approximately met
% ESDmin   = minimum modelled size
% ESDmax   = maximum modelled size
% ESDP1    = cell diameter of smallest prey class
% deltaOpt = optimal pred:prey diameter ratio for grazing

if ESDP1 > ESDmin
    
    n = double(n);
    
%     clear ESD_
    ESD_(1) = ESDP1; % create baseline grid -- lowest resolution compatible with ESDP1 and deltaOpt
    ESD_(2:n) = ESDP1 .* deltaOpt .^ (1:n-1);
    ESD_ = ESD_(ESD_ <= ESDmax); % discard sizes outside modelled range
    
    nc = length(ESD_) - 1; % number of 'complete' intervals within basline grid
    es = log10(ESDmax) - log10(ESD_(end)); % extra space beyond largest baseline-grid size
    r = round((n - 1) / (nc + es / log10(deltaOpt))); % grid resolution -- number of size classes per *deltaOpt diameter gap
    inc = log10(deltaOpt) / r; % diameter increments
    np = r * nc + 1 + floor(es / inc); % number of size classes -- best match to 'np' input
    
%     clear ESDP
    ESDP = nan(np, 1); % prey diameters
    ESDP(1) = ESDP1;
    ESDP(2:np) = ESDP1 .* deltaOpt .^ ((1:np-1) / r);
    
%     clear ESDZ
    
    ESDZ = ESDP; % predator diameters
    % When the predator cell diameters match the prey, the smallest
    % predators cannot optimally graze on any prey class -- unless the
    % prey preference function is normalised, but the whole purpose of this 
    % size class set up is to avoid normalising the preference function.
    
    % Uncomment to remove smallest predator sizes
%     ESDZ = ESDP(find(ESDP == ESD_(2)):end);
    % This choice ensures that all predators can optimally graze a prey class,
    % and that all size classes are within [ESDmin, ESDmax]. However, not all
    % prey classes are optimally grazed because predators are not large enough.
    
    % Larger predator size classes could be included to ensure all prey
    % size classes are optimally grazed by some predator. However, these larger
    % size classes would exceed the ESDmax value, which is based on observation
    % -- and if predator classes are also subject to grazing then the problem
    % just extends because the largest predators would not be optimally grazed...
    % leave this line commented out.
    
    % ESDZ = [ESDZ, 10 .^ (log10(ESDP(end)) + (1:r) .* (log10(deltaOpt) / r))];
    
    VolP = pi / 6 * (ESDP) .^ 3;
    VolZ = pi / 6 * (ESDZ) .^ 3;
    
    % interval edges    
    ESDPedges = intervalEdges(ESDP, ESDmin, ESDmax);
    ESDZedges = intervalEdges(ESDZ, ESDmin, ESDmax);

    
    
else
    error('Diameter of smallest cells, ESDP1, must exceed modelled minimum, ESDmin.')
end

end

%%

function out = intervalEdges(x, xmin, xmax)
nx = length(x);
lx = log10(x);
inc = diff(lx(1:2));
if isrow(x) || iscolumn(x)
    out = nan(1, nx + 1);
    if iscolumn(x), out = out'; end
    
    out(2:end-1) = 0.5 .* (lx(1:end-1) + lx(2:end));
    out(1) = max(log10(xmin), out(2) - inc);
    out(end) = min(log10(xmax), out(end-1) + inc);
    out = 10 .^ out;
else
    out = nan;
    warning('intervalEdges.m input "x" must be a row or column vector')
end
end
