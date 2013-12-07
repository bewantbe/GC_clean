% Multi-variate time series Granger causality test in time domain
% this version use '\' operator (QR decomposition) to perform LSE
% much stable than pos_nGrangerT and pos_nGrangerT2
% but cost much more memory and more time to calculate

function [GC, Deps, Aall] = pos_nGrangerT_qr(X, m)
if (nargin ~= 2 && nargin ~= 3)
    error('Usage: [GC, Deps] = nGrangerT3(X, m)');
end
if (exist('showcond', 'var')==0)
    showcond = 0;
end

[p, len] = size(X);
% make sure that the average of X is zero
X = bsxfun(@minus, X, mean(X,2));

Z = zeros(m*p, len-m);
for l = 1 : m
    Z(1+(l-1)*p:l*p, :) = X(:, m+1-l:len-l);
end

Aall = -(X(:, m+1:len) / Z);

fcov0 = @(x) x*x'/(len-m);
Deps = fcov0(X(:, m+1:len) + Aall*Z);

if (p == 1)
    GC = 0;
    return ;
end

Depsj = zeros(p-1, p);        % echo column corresponding to a excluded variate
idx = true(m*p, 1);           % the index of lines we want to solve
idx(1:p:m*p) = false;
for k = 1 : p
    Acj = -(X(idx(1:p), m+1:len) / Z(idx, :));     % solve p-1 variates regression
    Depsj(:, k) = diag(fcov0(X(idx(1:p), m+1:len) + Acj*Z(idx, :)));
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
