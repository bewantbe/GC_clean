% Calculate circlar covariance matrixes
% call in the form
%   R = getcov_circle(X, m)
%   calculate
%    R(k) = X(t) * X(t-k)'
% call in the form
%   R = getcov_circle(X, Y, m)
%   calculate
%   R(k) = X(t) * Y(t-k)'

% time cost is about: O(len * m * p^2)

function R = getcov_circle(X, Y, m)
if (~exist('m','var') && exist('Y','var')) && length(Y)==1
  m = Y;
end
if (exist('X','var')==0 || exist('m','var')==0 || isempty(X)==1)
    disp('usage: R = getcov_circle(X, m)');
    error('X must be a p*len dim non-empty matrix');
end
[p, len] = size(X);
if (len-m < 1)                    % max possible value of m is len-1
    error('Function getcov: order m is too high!');
end
X = bsxfun(@minus, X, mean(X,2));
R = zeros(p, p*(m+1));
if length(Y)==1
  for k = 0 : m
      R(:, 1+p*k : p+p*k) = [X(:, k+1:len), X(:, 1:k)] * X';
  end
else
  for k = 0 : m;
      R(:, 1+p*k : p+p*k) = [X(:, k+1:len), X(:, 1:k)] * Y';
  end
end
R = R / len;
end
