function [alpha, mu, c] = fitDirichlet(x)
% Input a (K,N) matrix, x, of K-simplexes. Each of the j=1,...,N simplexes
% is assumed to be distributed: x_j ~ Dir(alpha). Thus x represents a
% random sample of N, Dir(alpha) distributed, simplexes.
% Function returns estimate of vector concentration parameter, alpha.
K = size(x, 1); % number of categories
N = size(x, 2); % number of samples
mu = mean(x, 2); % expectation parameter
mu = mu ./ sum(mu);
% lmu = log(mu);
log_x = log(x);
ltp = mean(log_x, 2);
% c = -0.5 * (K - 1) / sum(mu .* (ltp - lmu)); % initial guess for concentration (precision) scalar parameter
c = K; % initial guess for concentration (precision) scalar parameter

% algorithm to maximise log-likelihood wrt c
stop = false;
maxiter = 1e3;
tol = 1e-2; % require at least 100*tol% change in value to continue
j = 0;
ch = nan(maxiter, 1);
while ~stop
    % This works fairly well... struggles with very low c
    j = j + 1;
    ch(j) = c;
    g = N * (psi(0, c) - sum(mu .* psi(0, c * mu)) + sum(mu .* ltp)); % gradient of log-likelihood wrt c, holding mu fixed
%     gt(j) = g;
    gg = N * (psi(1, c) - sum(mu .^ 2 .* psi(1, c * mu))); % 2nd order gradient
    u = abs(1 / c + 1 / c ^ 2 * g / gg);
    cnew = 1 / u;
    c = cnew;
    if isnan(c), break; end
    if g == 0
        stop = true;
    end
    if j > 1
        z = abs(ch(j) - ch(j-1)) / ch(j-1);
        if z < tol
            stop = true;
        end
    end
    if j >= maxiter
        stop = true;
    end
end

alpha = c * mu;
