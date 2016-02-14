% Multi-variate conditional Granger causality calculation in time domain.  v1.2
% Use only covariance series to calculate GC.
% Covariance series R is a p * (p*(m+1)) dimension matrix.
%
% Time cost is about: O(p^4 * m^3)
% RAM  cost is about: O(2*(p*m)^2) * 8Byte

function [GC, Deps, Aall] = RGrangerT(R)

[p, m] = size(R);
m = round(m/p)-1;                 % fitting order

% covariance series that has negative time
RR = zeros(p,p*(2*m-1));
RR(:, m*p-p+1:end) = R(:,1:m*p);
for k = 1 : m-1
    RR(:,(m-k)*p-p+1:(m-k)*p) = R(:,k*p+1:k*p+p)';
end
% construct the big covariance matrix
covz = zeros(p*m, p*m);
for k = 0 : m-1
    covz(k*p+1:k*p+p, :) = RR(:, (m-k)*p-p+1 : (m+m-k)*p-p );
end

Rb = -R(:,1+p:p+p*m)';            % nonhomogeneous item
Aall = (covz \ Rb)';              % solve all-jointed regression, covz * Aall' = Rb
Deps = R(:,1:p) - Aall*Rb;        % variance matrix of noise term

if (sum(Deps(eye(p)==1)<0) > 0)
    warning('nGranger:BAD_COV','Bad condition problem occured, the answer is wrong. consider use pos_nGeangerT2');
    Deps(Deps<0 & eye(p)==1) = 5e-324;     % make sure diag are positive (TODO: Is there better idea to fix this?)
end

if (p == 1)
    GC = 0;
    return ;
end

Depsj = zeros(p-1, p);            % echo column corresponding to a excluded variate
idx0 = 1:p:m*p;                   % indexes of variate to exclude
for k = 1 : p
    idx = true(1, m*p);           % the index of lines we want to solve
    idx(idx0) = false;
    idxb= idx(1:p);
    Rbj = Rb(idx, idxb);
    Acj = (covz(idx, idx) \ Rbj)';     % solve p-1 variates regression
    Depsj(:, k) = diag(R(idxb, idxb) - Acj*Rbj);
    idx0 = idx0 + 1;              % next variate
end

if (sum(Depsj<0) > 0)
    warning('nGranger:BAD_COV','Bad condition problem occured, the answer is wrong. consider use pos_nGeangerT2');
    Depsj(Depsj<0) = 5e-324;      % make sure diag are positive (TODO: Is there better idea to fix this?)
end

% use variances to calculate GC
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
