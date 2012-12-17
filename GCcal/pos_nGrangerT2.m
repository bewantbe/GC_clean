% Multi-variate time series Granger causality test in time domain  v1.2
% this version use completely symmetric positive definite covx from Xt
% time and memory cost is little bigger than nGrangerT.m

% time   cost is about: O(len * p * m + p^3 * m^2)
% memory cost is about: O(len * p + (m*p)^2)

function [GC, Deps, Aall] = pos_nGrangerT2(X, m, showcond)
if (nargin ~= 2 && nargin ~= 3)
    error('Usage: [GC, Deps] = nGrangerT2(X, m [, showcond])');
end
if (exist('showcond', 'var')==0)
    showcond = 0;
end

[p, len] = size(X);
% make sure that the average of X is zero
%X = X - mean(X,2)*ones(1,len);
X = bsxfun(@minus, X, mean(X,2));

t_bg = m+1;     % time range in common
t_ed = len-m;
% middle part
R = zeros(p, (2*m+1)*p);
for k = 0 : m
    R_k = X(:, t_bg:t_ed) * X(:, t_bg-k:t_ed-k)';
    R(:, (m+k)*p+1 : (m+1+k)*p) = R_k;
    R(:, (m-k)*p+1 : (m+1-k)*p) = R_k';
end
% head and tail parts
covz = zeros((m+1)*p,(m+1)*p);
for i1=0:m
    for i2=0:m
        k=i2-i1;
        if (k>0)
            covz(i1*p+1:(i1+1)*p, i2*p+1:(i2+1)*p) =...
              X(:,m+1-i1:t_bg-1)*X(:,m+1-i1-k:t_bg-1-k)'...
            + X(:,t_ed+1:len-i1)*X(:,t_ed+1-k:len-i1-k)';
        else
            covz(i1*p+1:(i1+1)*p, i2*p+1:(i2+1)*p) =...
              X(:,m+1-i2+k:t_bg-1+k)*X(:,m+1-i2:t_bg-1)'...
            + X(:,t_ed+1+k:len-i2+k)*X(:,t_ed+1:len-i2)';
        end
    end
end
% combine them
for k=0:m
    covz(k*p+1:(k+1)*p,:) = covz(k*p+1:(k+1)*p,:) + R(:,(m-k)*p+1:(2*m+1-k)*p);
end
covz = covz/(len-m);
R = zeros(p, (m+1)*p);
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

if (p == 1)
    GC = 0;
    return ;
end

Depsj = zeros(p-1, p);           % echo column corresponding to a excluded variate
idx0 = 1:p:m*p;                  % indexes of variate to exclude
for k = 1 : p
    idx = true(1, m*p);          % the index of lines we want to solve
    idx(idx0) = false;
    idxb= idx(1:p);
    Rbj = Rb(idx, idxb);
    Acj = (covz(idx, idx) \ Rbj)';     % solve p-1 variates regression
    Depsj(:, k) = diag(R(idxb, idxb) - Acj*Rbj);
    idx0 = idx0 + 1;             % next variate
end

%Depsj(Depsj<0) = 0; % make sure diag are positive (TODO: Is there better idea to fix this?)

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
