% Give order choice for autoregression and joint-regression from covariance

function [od_joint, od_vec] = chooseROrderFull(R, len, ic_mode, od_max)
if nargin < 2
  error('Usage: [od_joint, od_vec] = chooseROrderFull(R, len, ic_mode, od_max)');
end
if ~exist('ic_mode', 'var')
  ic_mode = 'AIC';
end
slow_mode = 1;

p = size(R,1);
od_max = round(size(R,2)/p)-1;

od_joint = chooseROrder(R, len, ic_mode, slow_mode);
od_vec = zeros(1,p);
for id_y = 1:p
  idx0 = true(1,p);  idx0(id_y)=false;
  idx  = true(1,(od_max+1)*p);
  idx(idx0(rem((1:((od_max+1)*p))-1,p)+1)==false) = false;   % exclude all id_y
  od_auto = chooseROrder(R(idx0,idx), len, ic_mode, slow_mode);
  od_vec(id_y) = od_auto;
end

end
