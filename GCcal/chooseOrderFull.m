% Give order choice for autoregression and joint-regression

function [od_joint, od_vec] = chooseOrderFull(X, ic_mode, od_max)
if nargin < 1
  error('Usage: [od_joint, od_vec] = chooseOrderFull(X, ic_mode, od_max)');
end
if ~exist('ic_mode', 'var')
  ic_mode = 'AIC';
end
if ~exist('od_max', 'var')
  od_max = 499;
end

[p, len] = size(X);
R = getcovpd(X, od_max);

[od_joint, od_vec] = chooseROrderFull(R, len, ic_mode, od_max);

end
