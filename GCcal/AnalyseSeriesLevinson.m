% Analyse a given time series in a series of orders
% Input & Output: See comments in AnalyseSeries.m
% Since this function relies on Levinson algorithm, it not suitable for
% bad condition problem (e.g. very short data, filtered data)
%
% Time cost:
% O(p^2*m*L) + O(p^3*m^2) + O( m*(p^2*m^2*log(m)) ) + O( m*(p*m^3 + p^2*m^2) )
%     cov       Levinson       (partial)inversion                 GC

function [oGC, oDe, R] = AnalyseSeriesLevinson(X, s_od, b_input_cov)
  if exist('b_input_cov','var') && b_input_cov
    % so X is actually covariance
    R = X;
    p = size(R, 1);
    m = round(size(R, 2) / p) - 1;
    if (m < max(s_od))
      error('GC_clean: AnalyseSeriesLevinson: no enough data');
    end
  else
    p = size(X, 1);
    m = max(s_od);
    R = getcovpd(X, m);
  end
  oGC = zeros(p, p, length(s_od));
  oDe = zeros(p, p, length(s_od));

  % Prepare for Levinson algo.
  % Reverse order R: [R(-m),...,R(-1),R(0)]. Note R(-k)=R(k)'
  b_R = reshape(flipdim(...
          permute(reshape(R,p,p,[]),[2,1,3])...
          , 3), p, []);
  R3d = reshape(R,p,p,[]);
  f_row2col = @(x) reshape(permute(reshape(x,p,p,[]), [2,1,3]),p,[]).';
  % RHS: covz * xf = yf, covz * xb = yb
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
  if sum(s_od==1) > 0  % we need GC at this order
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
    % Got solution xf, xb for
    %   T(1:p*k, 1:p*k) * xf = yf(1:p*k, :)
    %   T(1:p*k, 1:p*k) * xb = yb(end-p*k+1:end, :)
    if sum(s_od==k) > 0
      kk = kk + 1;
      [oGC(:,:,kk), oDe(:,:,kk)] = GCfromQyyA(R, fk, bk, xf, xb);
    end
  end
end

% see also RGrangerTfast.m
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
  A = permute(reshape(-A, p,p,[]), [3,1,2]);   % index: (time_lag, i, j)  (j->i)
  d = diag(D);
  GC = zeros(p, p);
  for j = 1 : p
    GC(:, j) = log1p(sum(QyyArr(:,:,j) \ A(:,:,j) .* A(:,:,j),1)' ./ d);
    GC(j, j) = 0;
  end
end

% Get diagonal block component of inverse block Levinson matrix
% Expample
%   [cf3d, cb3d, df3d, db3d, A] = BlockLevinson_v3_inv(R, 1);
%   QyyArr = GetDiagInvCov(cf3d, cb3d, df3d, db3d)
% Equivalent to
%   covz = R2covz(R);
%   covz1 = covz(p+1:end, p+1:end);
%   Qz = inv(covz1);
%   id_0 = 0:p:p*m-1;
%   for j = 1 : p
%     QyyArr(:,:,j) = Qz(id_0+j, id_0+j);
%   end

function QyyArr = GetDiagInvCov(cf3d, cb3d, df3d, db3d)
  p = size(cf3d,1);
  m = size(cf3d,3);

  if m==1
    QyyArr = reshape(diag(cf3d),1,1,[]);
    return;
  end

  cf3d = permute(cf3d, [2,1,3]);  % use the transposed form
  df3d = permute(df3d, [2,1,3]);  % use the transposed form

  QyyArr = zeros(p,m,m);
  f_T2 = fft(db3d(:,:,[end, 1:end-1]), m, 3);
  f_T4 = fft(cb3d(:,:,[end, 1:end-1]), m, 3);
  df3d(:,:,1) = df3d(:,:,1) - eye(p);
  % The FFT in Octave(3.6.2) can not do fft along singleton 3rd dimension.
  % Like fft(rand(2,2,1), 9, 3)
  j = 1;
  f_T1 = fft(cat(3,cf3d(:,:,1),zeros(p,p,1)), m, 3);
  f_T3 = fft(cat(3,df3d(:,:,1),zeros(p,p,1)), m, 3);
  QyyArr(:,:,j) = squeeze(...
    sum(f_T1.*f_T2, 1) -...
    sum(f_T3.*f_T4, 1));
  QyyArr(:,:,j) = real(ifft(QyyArr(:,:,j), m, 2));

  for j=2:m
    QyyArr(:,:,j) = real(ifft(squeeze(...
      sum(fft(cf3d(:,:,j:-1:1), m, 3).*f_T2) -...
      sum(fft(df3d(:,:,j:-1:1), m, 3).*f_T4)...
      ), m, 2));
  end
  QyyArr = permute(QyyArr, [3,2,1]);
end

