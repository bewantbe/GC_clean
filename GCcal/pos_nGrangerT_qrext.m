% Multi-variate conditional Granger causality calculation in time domain.
% This version use '\' operator (QR decomposition) to perform LSE.
% Much stable than pos_nGrangerT and pos_nGrangerT2,
% but cost much more memory and more time to calculate.
% This version allow user put time series `extX' to be project out.
%
% Time cost: O(p * (p*m+extl)^2 * L),  extl is number of variables of extX
% RAM  cost: O((p*m+extl)*L)

function [GC, Deps, Aall, res] = pos_nGrangerT_qrext(X, m, extX)
if (nargin ~= 3)
    error('Usage: [GC, Deps] = pos_nGrangerT_qrext(X, m, extX)');
end

[p,   len] = size(X);
[p2, len2] = size(extX);
if (len ~= len2)
  error('length of extX and X should be the same');
end

Z = zeros(m*p+p2, len-m);  % the last row is for extX
for l = 1 : m
    Z(1+(l-1)*p:l*p, :) = X(:, m+1-l:len-l);
end
Z(m*p+1:end, :) = extX(:, m+1:len);

Aall = -(X(:, m+1:len) / Z);

fcov0 = @(x) x*x'/(len-m-1);
Deps = fcov0(X(:, m+1:len) + Aall*Z);

res = X(:, m+1:len) + Aall*Z;

if (p == 1)
    GC = 0;
    return ;
end

Depsj = zeros(p-1, p);     % echo column corresponding to a excluded variate
idx = true(m*p, 1);        % the index of lines we want to solve
idx(1:p:m*p) = false;
for k = 1 : p
    % solve p-1 variable regression (without k-th variable)
    Acj = -(X(idx(1:p), m+1:len) / Z([idx;true(p2,1)], :));
    Depsj(:, k) = diag(fcov0(X(idx(1:p), m+1:len) + Acj*Z([idx;true(p2,1)], :)));
    idx = circshift(idx,1);
end

GC = zeros(p,p);
dd = diag(Deps);
dj = zeros(p,1);
for k = 1 : p
    dj(1:p-1) = Depsj(:, k);
    dj(k+1:end) = dj(k:p-1);
    dj(k) = dd(k);
    GC(:, k) = log(dj ./ dd);
end

end
