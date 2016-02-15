% Multi-variate conditional Granger causality calculation in time domain.
% Much faster than RGrangerT(R) and pos_nGrangerT2() for large (>50) variable case.
% Almost as stable as pos_nGrangerT2(), and mathematically equivalent to it.
% Time cost: O( (p*m)^3 )              (Typical time for p=100,m=40 is 5 seconds)
% RAM  cost: O( 3.5*(p*m)^2 ) * 8Byte  (Typical max scale is p=400,m=40 for 8GB RAM)

function [GC, D, A2d] = RGrangerTfast(R, p)
  if exist('p', 'var')
    % so R is covz matrix
    covz = R;
    [A2d, D] = ARregressionpd(covz, p);
  else
    if size(R, 1) == size(R, 2)
      error('Dimension incompatible. Do you mean RGrangerTfast(covz, p)');
    end
    % so R is 2d cov series
    p = size(R, 1);
    covz = R2covz(R);
    [A2d, D] = ARregressionpd(covz, p);
  end
  m = size(covz, 2) / p - 1;

  Qz = inv(covz(p+1:end, p+1:end));
  a = reshape(-A2d, p,p,[]);
  a = permute(a, [3,1,2]);   % index: (time_lag, i, j)  (j->i)
  d = diag(D);
  id_0 = 0:p:p*m-1;

  GC = zeros(p, p);
  for j = 1 : p
    Qjj = Qz(id_0+j, id_0+j);
    GC(:, j) = log1p(sum(Qjj \ a(:,:,j) .* a(:,:,j))' ./ d);
    GC(j, j) = 0;
  end
end
