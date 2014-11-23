% Analysis a given time series in different order
% Input & Output: See comments in AnalyseSeries.m
% Since this function relies on Levinson algorithm, it not suitable for
% bad condition problem (e.g. very short data, filtered data)

% Time cost (haven't verified):
% O(p^2*m*L) + O(p^3*m^2*log(m))

function [oGC, oDe, R] = AnalyseSeriesLevinson(X, s_od)
  if ~exist('bad_mode','var')
    bad_mode = 0;
  end
  p = size(X, 1);
  m = max(s_od);
  R = getcovpd(X, m);
  oGC = zeros(p, p, length(s_od));
  oDe = zeros(p, p, length(s_od));

  % Prepare for Levinson algo.
  % The 2d form [R(-m),...,R(-1),R(0)]. Note R(-k)=R(k)'
  b_R = reshape(flipdim(...
          permute(reshape(R,p,p,[]),[2,1,3])...
          , 3), p, []);
  R3d = reshape(R,p,p,[]);
  f_row2col = @(x) reshape(permute(reshape(x,p,p,[]), [2,1,3]),p,[]).';
  % RHS
  yf = reshape(R3d(:,:,2:end),p,[])';
  yf = [yf, [zeros(p,p); yf(1:p*(m-1),:)]];
  yb = reshape(permute(R3d(:,:,end:-1:2),[2,1,3]),p,[]).';
  yb = [yb, [yb(p+1:p*m,:); zeros(p,p)]];

  % first Levinson step
  T1 = R(:, 1:p);
  fk = inv(T1);
  bk = fk;
  xf = T1 \ yf(1:p, :);
  xb = T1 \ yb(end-p+1:end, :);
  p2f = size(yf,2);
  p2b = size(yb,2);

  kk = 0;
  if sum(s_od==1) > 0
    kk = kk + 1;
    [oGC(:,:,kk), oDe(:,:,kk)] = GCfromQyyA(R, fk, bk, xf, xb);
  end
  for k=2:m
    % Levinson steps
    eps_fk = b_R(:, end-k*p+1:end-p) * fk;
    eps_bk = R(:, p+1:k*p) * bk;
    Vfb = [[fk; zeros(p)], [zeros(p); bk]] / [eye(p), eps_bk; eps_fk, eye(p)];
    fk = Vfb(:,1:p);
    bk = Vfb(:,p+1:end);
    eps_xf = b_R(:, end-k*p+1:end-p) * xf(1:p*(k-1), :);
    xf = [xf; zeros(p,p2f)] + bk * (yf(p*(k-1)+1:p*k, :) - eps_xf);
    eps_xb = R(:, p+1:k*p) * xb(1:p*(k-1), :);
    xb = [zeros(p,p2b); xb] + fk * (yb(end-p*k+1:end-p*(k-1), :) - eps_xb);
    if sum(s_od==k) > 0
      kk = kk + 1;
      [oGC(:,:,kk), oDe(:,:,kk)] = GCfromQyyA(R, fk, bk, xf, xb);
    end
    %norm(T(1:p*k, 1:p*k) * xf - yf(1:p*k, :) )
    %norm(T(1:p*k, 1:p*k) * xb - yb(end-p*k+1:end, :) )
  end

end

function [GC, D] = GCfromQyyA(R, fk, bk, xf, xb)
  p = size(fk, 2);
  m = size(fk, 1) / p;

  f_2Dto3D  = @(x) permute(reshape(x,p,m,[]),[1,3,2]);
  f_2Dto3DT = @(x) permute(reshape(x,p,m,[]),[3,1,2]);
  fk = f_2Dto3D (fk);
  bk = f_2Dto3DT(bk);
  df = f_2Dto3D (xf(:, p+1:p+p) - xb(:, 1:p));
  db = f_2Dto3DT(xb(:, p+1:p+p) - xf(:, 1:p));
  QyyArr = GetDiagInvCov(fk, bk, df, db);
  A = -xf(:, 1:p)';
  D = R(:,1:p*(m+1)) * [eye(p,p); A'];         % the residual variance
  A = permute(reshape(-A, p,p,[]), [3,1,2]);
  d = diag(D);
  GC = zeros(p, p);
  for j = 1 : p
    GC(:, j) = log1p(sum(QyyArr(:,:,j) \ A(:,:,j) .* A(:,:,j))' ./ d);
    GC(j, j) = 0;
  end
end
