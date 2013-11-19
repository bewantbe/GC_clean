% Generate a p dimension multivariate time series

% Call this function with following parameters will
% give a pre-define AR model:
%   gdata(len, 2, 2)
%   gdata(len, 2, 3)
%   gdata(len, 3, 5)

% len: length of series
% m  : Order of model
% p  : number of variates
% X(t) is a p dimension column vector, t=1..len
% noisecov is the covariance matrix of noise terms(must symmetric positive definite)
% A is the coefficient matrices, p*(p*m) dimension

% model: X(t) + A1*X(t-1) + A2*X(t-2) +...+ Am*X(t-m) = noise(t)

% time cost is about: O(len * m * p^2)

function X = gdata(len, m, p, noisecov, A)
%tic();
if (exist('len', 'var') == 0)
    len = 100;
end
if (exist('m','var') == 0)
    m = 2;
end
if (exist('p','var') == 0)
    p = 2;
end
if (exist('noisecov','var') == 0)
    if (p==2)
        noisecov = [1.0, 0.4; 0.4, 0.7];
    end
    if (p==3)
        noisecov = diag([0.3 1.0 0.2]);
    end
    if (p==5)
        noisecov = diag([0.6 0.5 0.3 0.3 0.6]);
    end
end
if (exist('A','var') == 0)
    if (m==2 && p==2)
        A = [-0.9 ,  0.0, 0.5, 0.0;
             -0.16, -0.8, 0.2, 0.5];
    end
    if (m==2 && p==3)
        A = [-0.8  0.0 -0.4  0.5 -0.2  0.0;
              0.0 -0.9  0.0  0.0  0.8  0.0;
              0.0 -0.5 -0.5  0.0  0.0  0.2];
    end
    if (m==3 && p==5)
        s2 = sqrt(2);
        A = zeros(5, 5*3);
        A(1,1) =-0.95*s2;  A(1,6) = 0.9025;
        A(2,6) =-0.5;
        A(3,11)= 0.4;
        A(4,6) = 0.5;      A(4,4) =-0.25*s2;  A(4,5) =-0.25*s2;
        A(5,4) = 0.25*s2;  A(5,5) =-0.25*s2;
    end
end
extnum = 200;                           % wait for stable state
len = len + extnum;
% generate noise terms
tf = chol(noisecov);                    % Compute the Cholesky factor
Eps = tf' * randn(p, len);              % Noise terms
X = zeros(p, len);                      % pre allocate memory
for t = 1 : len
    mmax = min(t-1, m);             % avoid X(t) out of range
    X(:,t) = Eps(:,t);
    for k = 1 : mmax
        X(:,t) = X(:,t) - A(:, 1+p*(k-1):p+p*(k-1)) * X(:,t-k);
    end
    %X(:,t) = Eps(:,t);
    %X(:,t)-= X()*;                  % $!#%$%$&%^*&^>(
end
X(:,1:extnum) = [];
%disp('cov matrix of noise term');
%disp(cov(Eps'));    % not assume ave = 0
%disp(Eps*Eps'/len);   % assume ave = 0
%disp('in gdata');
%toc();
end
