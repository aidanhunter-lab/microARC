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
    Vol_low_mu_max_lim = FixedParams.PPsize(1);
end
if ~exist('Vol_high_mu_max_lim', 'var')
    Vol_high_mu_max_lim = FixedParams.PPsize(end);
end
if ~exist('mu_max_lowest', 'var')
    mu_max_lowest = 0.2;
end

i_Qmin_a = strcmp('Qmin_QC_a', parNames);
i_Qmin_b = strcmp('Qmin_QC_b', parNames);
i_Vmax_a = strcmp('Vmax_QC_a', parNames);
i_Vmax_b = strcmp('Vmax_QC_b', parNames);
i_pmax_a = strcmp('pmax_a', parNames);
i_pmax_b = strcmp('pmax_b', parNames);

Params.constraints = @(x) all_constraints(x, i_Qmin_a, i_Qmin_b, ...
    i_Vmax_a, i_Vmax_b, i_pmax_a, i_pmax_b, ...
    Vol_low_mu_max_opt, Vol_high_mu_max_opt, ...
    Vol_low_mu_max_lim, Vol_high_mu_max_lim, mu_max_lowest);

end

%--------------------------------------------------------------------------

function v = vol_mu_max_opt(x, i_Qmin_a, i_Qmin_b, i_Vmax_a, i_Vmax_b, ...
    i_pmax_a, i_pmax_b)
% Return the log volume at optimal growth rate, mu_max.
v = (log(x(i_pmax_a)) + log(x(i_Qmin_a)) + log(x(i_Vmax_b) - x(i_Qmin_b)) - ...
    log(-x(i_pmax_b)) - log(x(i_Vmax_a))) / ...
    (x(i_Vmax_b) - x(i_pmax_b) - x(i_Qmin_b));
end

function c = constraint_mu_max_opt(x, iqa, iqb, iva, ivb, ipa, ipb, vol_low, vol_high)
vol_opt = vol_mu_max_opt(x, iqa, iqb, iva, ivb, ipa, ipb);
 c1 = log(vol_low) - vol_opt;
 c2 = vol_opt - log(vol_high);
 c = [c1, c2];
end

function c = constraint_mu_max_lim(x, iqa, iqb, iva, ivb, ipa, ipb, vol_low, vol_high, mu_max_lb)
v = [vol_low, vol_high];
m = x(ipa) .* x(iva) .* (v .^ (x(ipb) + x(ivb))) ./ ...
    (x(ipa) .* x(iqa) .* (v .^ (x(ipb) + x(iqb))) + x(iva) .* (v .^ x(ivb)));
c = mu_max_lb - m;
end

function [c, ceq] = all_constraints(x, iqa, iqb, iva, ivb, ipa, ipb, ...
    vol_low_opt, vol_high_opt, vol_low_lim, vol_high_lim, mu_max_lb)
c_opt = constraint_mu_max_opt(x, iqa, iqb, iva, ivb, ipa, ipb, vol_low_opt, vol_high_opt);
c_lim = constraint_mu_max_lim(x, iqa, iqb, iva, ivb, ipa, ipb, vol_low_lim, vol_high_lim, mu_max_lb);
c = [c_opt(:); c_lim(:)];
ceq = [];
end

