% Multi-variate conditional Granger causality calculation in time domain.
% This version use '\' operator (QR decomposition) to perform LSE,
% Much stable than pos_nGrangerT and pos_nGrangerT2,
% but cost much more memory and more time to calculate.
%
% Time cost: O(p^3*m^2*L)
% RAM  cost: O(p*m*L)

function [GC, Deps, Aall] = pos_nGrangerT_qrm(X, m)
if (nargin ~= 2 && nargin ~= 3)
    error('Usage: [GC, Deps] = nGrangerT3(X, m)');
end
if (exist('showcond', 'var')==0)
    showcond = 0;
end

[p, len] = size(X);
% make sure that the average of X is zero
%X = bsxfun(@minus, X, mean(X,2));

Z = zeros(m*p+1, len-m);  % the last row is for mean value
for l = 1 : m
    Z(1+(l-1)*p:l*p, :) = X(:, m+1-l:len-l);
end
Z(m*p+1, :) = 1;

Aall = -(X(:, m+1:len) / Z);

fcov0 = @(x) x*x'/(len-m-1);
Deps = fcov0(X(:, m+1:len) + Aall*Z);

if (p == 1)
    GC = 0;
    return ;
end

Depsj = zeros(p-1, p);        % echo column corresponding to a excluded variate
idx = true(m*p, 1);           % the index of lines we want to solve
idx(1:p:m*p) = false;
for k = 1 : p
    % solve p-1 variable regression (without k-th variable)
    Acj = -(X(idx(1:p), m+1:len) / Z([idx;true], :));
    Depsj(:, k) = diag(fcov0(X(idx(1:p), m+1:len) + Acj*Z([idx;true], :)));
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
