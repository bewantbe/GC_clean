% Calculate high order regression of AR model
% by using covariance R

function [Aall, Deps] = ARregressionpd(R, p)
if exist('p','var') && size(R,1)==size(R,2)
  m = size(R,1)/p-1;
  covz = R(p+1:end, p+1:end);
  R = R(1:p, 1:end);
else
  if exist('p','var')
    warning('ARregressionpd(R, p): input parameter "p" ignored');
  end
  [p, m] = size(R);
  m = round(m/p)-1;

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
end

Rb = -R(:,1+p:p+p*m)';           % nonhomogeneous item
Aall = (covz \ Rb)';             % solve all-jointed regression, covz * Aall' = Rb
Deps = R(:,1:p) - Aall*Rb;       % variance matrix of noise term
if (p>1)
  Deps(Deps<0 & eye(p)==1) = 0;  % make sure their are positive (TODO: Is there better idea to fix this?)
end

end
