% Multi-variate time series Granger causality test in time domain  v1.1

% X is p*len dim matrix, means p-variate data length len
% m is the model order desired to estimate
% GC is the causality matrix, influence direction is column -> row
% Deps is the variance matrix of noise term (D(epsilon))

% time cost is about: O(len * m * p^2 + p^3 * m^2)

function [GC, Deps, Aall] = nGrangerT(X, m)
if (nargin ~= 2)
    error('Usage: nGrangerT(X, m), m is the order of AR');
end

[GC, Deps, Aall] = RGrangerT(getcovpd(X, m));

end
