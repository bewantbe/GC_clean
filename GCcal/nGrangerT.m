% Multi-variate conditional Granger causality calculation in time domain.  v1.1
% Var  Type    Size     Meaning
% X    matrix  p*len    means p-variate data length len
% m    scalar  1*1      model order desired to estimate
% GC   matrix  p*p      causality matrix, influence direction is column -> row
% Deps matrix  p*p      variance matrix of noise term (D(epsilon))
% Aall matrix  p*(p*m)  Fitting coef: [A(1), A(2), ..., A(m)]
%
% Time cost is about: O(len * m * p^2) + O(p^4 * m^3)

function [GC, Deps, Aall] = nGrangerT(X, m)
if (nargin ~= 2)
    error('Usage: nGrangerT(X, m), m is the order of AR');
end

[GC, Deps, Aall] = RGrangerT(getcovpd(X, m));

end
