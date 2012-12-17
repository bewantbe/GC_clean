% Compute pair-wise Granger Causality of p-variable input

% time cost is about: O(len * m * p^2 + p^2 *8*m^2)

function GC = pairGrangerT(X, m)
if (nargin ~= 2)
    error('Usage: pairGrangerT(X, m), m is the order of AR');
end

GC = pairRGrangerT(getcovpd(X, m));

end

