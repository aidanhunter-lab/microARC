function Params = parameterConstraints(FixedParams, Params, varargin)

% Non-linear constraints need to appear in argument after boundsUpper. This
% argument must be a function accepting X as an argument and returning
% vectors C and Ceq, where C(X) <= 0 and Ceq(X) = 0 define the constraints.
% Store these constraints in the Params struct.

extractVarargin(varargin)

parNames = FixedParams.tunePars;

% Constrain mu_max (see Ward et al. 2017) to have maximum for cell volume
% between Vol_low_mu_max_opt and Vol_high_mu_max_opt
if ~exist('Vol_low_mu_max_opt', 'var')
    Vol_low_mu_max_opt = 10;
end
if ~exist('Vol_high_mu_max_opt', 'var')
    Vol_high_mu_max_opt = 100;
end
if ~exist('Vol_low_mu_max_lim', 'var')
%     Vol_low_mu_max_lim = FixedParams.PPsize(1);
    Vol_low_mu_max_lim = 1e-1; % mu m^3
end
if ~exist('Vol_high_mu_max_lim', 'var')
%     Vol_high_mu_max_lim = FixedParams.PPsize(end);
    Vol_high_mu_max_lim = 1e5;
end
if ~exist('mu_max_lowest', 'var')
    % Lower bound on mu_max at volumes [Vol_low_mu_max_lim, Vol_high_mu_max_lim]
%     mu_max_lowest = 0.2;
    mu_max_lowest = [0.4, 0.25];
end

constrainedPars = {'Qmin_QC_a' 'Qmin_QC_b' 'Vmax_QC_a' 'Vmax_QC_b' 'pmax_a' 'pmax_b'};
nc = length(constrainedPars);
ps.parName = constrainedPars;
i = contains(parNames, constrainedPars); % index constrained params that vary
ip = parNames(i);
ps.fixed = ~contains(ps.parName, ip);
ps.index = nan(1, nc);
ps.fixedVals = nan(1, nc);
for j = 1:nc
    if ~ps.fixed(j)
        ps.index(j) = find(strcmp(ps.parName{j}, parNames));
    else
        ps.fixedVals(j) = Params.(ps.parName{j});
    end
end

ps.Vol_low_mu_max_opt = Vol_low_mu_max_opt;
ps.Vol_high_mu_max_opt = Vol_high_mu_max_opt;
ps.Vol_low_mu_max_lim = Vol_low_mu_max_lim;
ps.Vol_high_mu_max_lim = Vol_high_mu_max_lim;
ps.mu_max_lowest = mu_max_lowest;

Params.constraints = @(x) all_constraints(x, ps);

end

%--------------------------------------------------------------------------

function p = get_constrained_params(x, s)
n = s.parName;
N = length(n);
y = nan(1, N);
y(~s.fixed) = x(s.index(~s.fixed));
y(s.fixed) = s.fixedVals(s.fixed);
for i = 1:N
    p.(n{i}) = y(i);
end
end

function v = vol_mu_max_opt(p)
% Return the log volume at optimal growth rate, mu_max.
v = (log(p.pmax_a) + log(p.Qmin_QC_a) + log(p.Vmax_QC_b - p.Qmin_QC_b) - ...
    log(-p.pmax_b) - log(p.Vmax_QC_a)) / ...
    (p.Vmax_QC_b - p.pmax_b - p.Qmin_QC_b);
end

function c = constraint_mu_max_opt(p, vol_low, vol_high)
vol_opt = vol_mu_max_opt(p);
 c1 = log(vol_low) - vol_opt;
 c2 = vol_opt - log(vol_high);
 c = [c1, c2];
end

function c = constraint_mu_max_lim(p, vol_low, vol_high, mu_max_lb)
v = [vol_low, vol_high];
m = p.pmax_a .* p.Vmax_QC_a .* (v .^ (p.pmax_b + p.Vmax_QC_b)) ./ ...
    (p.pmax_a .* p.Qmin_QC_a .* (v .^ (p.pmax_b + p.Qmin_QC_b)) + p.Vmax_QC_a .* (v .^ p.Vmax_QC_b));
c = mu_max_lb - m;
end

function [c, ceq] = all_constraints(x, s)
p = get_constrained_params(x, s);
% pn = s.parName;
c_opt = constraint_mu_max_opt(p, s.Vol_low_mu_max_opt, s.Vol_high_mu_max_opt);
c_lim = constraint_mu_max_lim(p, s.Vol_low_mu_max_lim, s.Vol_high_mu_max_lim, s.mu_max_lowest);
c = [c_opt(:); c_lim(:)];
ceq = [];
end



