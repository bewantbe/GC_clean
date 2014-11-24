% Calculate GC without explicitly do p-1 variable regressions
% Much faster than RGrangerT(R) for large variable case.

function GC = RGrangerTfast(R, p)
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
