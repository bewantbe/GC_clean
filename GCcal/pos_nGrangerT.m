% Multi-variate conditional Granger causality calculation in time domain.
% This version (1.0) use completely symmetric positive definite covz from Xt.
% Much slower and cost much more memory than nGrangerT.m.
% Recommend to use pos_nGrangerT2() instead.
%
% Time cost is about: O(len * m^2 * p^2) + O(p^4 * m^3)

function [GC, Deps, Aall] = pos_nGrangerT(X, m, showcond)
if (nargin ~= 2 && nargin ~= 3)
    error('Usage: nGrangerT(X, m)');
end
if (exist('showcond', 'var')==0)
    showcond = 0;
end

[p, len] = size(X);
% make sure that the average of X is zero
%X = X - mean(X,2)*ones(1,len);
X = bsxfun(@minus, X, mean(X,2));
% completely symmetric positive definite version
Z = zeros((m+1)*p,len-m);
for l = 0 : m
    Z(1+l*p:p+l*p,:) = X(:,m-l+1:len-l);
end
covz = Z*Z'/(len-m);
R(:,1:p*(m+1)) = covz(1:p,1:p*(m+1));
%covz = covz(1:m*p, 1:m*p);
covz = covz(p+1:p+m*p, p+1:p+m*p);

if showcond>0
    disp(['det  = ' num2str(det(covz))]);
    disp(['cond = ' num2str(cond(covz))]);
% there is no isdefinite() in matlab
%disp(['positive definite? ' num2str(isdefinite(covz))]);
end

Rb = -R(:,1+p:p+p*m)';        % nonhomogeneous item
Aall = (covz \ Rb)';          % solve all-jointed regression, covz * Aall' = Rb
Deps = R(1:p,1:p) - Aall*Rb;  % variance matrix of noise term

%save('-ascii','Aall.txt','Aall');    % save coefficient
%save('-ascii','Deps.txt','Deps');    % save noise term

idx0 = 1:p:m*p;               % index of variates to exclude

% variance matrix of noise term exclude 1 variate one by one
Depsj = zeros(p-1, (p-1)*p);
for l = 0 : p-1
    idx = zeros(1, m*p);
    idx(idx0) = 1;
    idx = find(idx==0);       % the index of lines we want to solve
    idxb= idx(1:p-1);
    Rbj = Rb(idx, idxb);
    covzj = covz(idx, idx);
    Acj = (covzj \ Rbj)';     % solve p-1 variates regression
    Depsj(:,1+(p-1)*l:p-1+(p-1)*l) = R(idxb, idxb) - Acj*Rbj;
    idx0 = idx0 + 1;
end

% causality matrix of k(column) -> l(row)
GC = zeros(p,p);
for l = 1 : p
    for k = 1 : p
        if (k==l)
            continue;
        end
        if (k>l)
            ij = l;
        else
            ij = l-1;
        end
        GC(l,k) = log( Depsj(ij, (p-1)*(k-1)+ij) / Deps(l,l) );
    end
end
end
