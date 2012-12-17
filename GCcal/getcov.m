% Calculate covariance matrices  v1.1

% X is a multivariate time series p*len dimension
% p is the number of variates, len is the length of time
% m is the maximun offset of covariance to calculate
% result R is a p * (p*(m+1)) dimension matrix

% time cost is about: O(len * m * p^2)

function R = getcov(X, m)

if (exist('X','var')==0 || exist('m','var')==0 || isempty(X)==1)
	disp('usage: R = getcov(X, m)');
	error('X must be a p*len dim non-empty matrix');
end
[p, len] = size(X);
if (len-m < 1)                     % max possible value of m is len-1
	error('Function getcov: order m is too high!');
end
%X = X - mean(X, 2)*ones(1,len);   % make sure that the average of X is zero
X = bsxfun(@minus, X, mean(X,2));  % 3 times slower in Octave 3.2.4, But OK in 3.4.2
R = zeros(p, p*(m+1));
for k = 0 : m
	R(:, 1+p*k : p+p*k) = (X(:, k+1:len) * X(:, 1:len-k)') / (len-k);
end

end
