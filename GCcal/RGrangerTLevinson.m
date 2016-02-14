% Multi-variate conditional Granger causality calculation in time domain.
% Do regression by Levinson method (Toeplitz)
% Significantly faster than RGrangerTfast(R) for large variables case, and use less RAM.
% Time cost: O( p^3*m^2*log(m) )    (Typical time for p=100,m=40 is 1.7 seconds)
% RAM  cost: O( 24*p^2*m ) * 8Byte  (Typical max scale is p=1000,m=40 for 8GB RAM)

function GC = RGrangerTLevinson(R)
  [s_Qjj, A, D] = GetQyy(R);

  p = size(R, 1);
  A = reshape(-A, p,p,[]);  % index: (i, j, time_lag)
  A = permute(A, [3,1,2]);  % index: (time_lag, i, j)  (j->i)
  d = diag(D);

  GC = zeros(p, p);
  for j = 1 : p
    GC(:, j) = log1p(sum(s_Qjj(:,:,j) \ A(:,:,j) .* A(:,:,j))' ./ d);
    GC(j, j) = 0;
  end
end


% Copy from local repo test_blocklevinson
function [QyyArr, A, D] = GetQyy(R)
  [p, m] = size(R);
  m = round(m/p) - 1;

  %R = [R, zeros(p,p)];
  [cf3d, cb3d, df3d, db3d, A] = BlockLevinson_v3_inv(R, 1);
  D = R * [eye(p,p); A'];         % the residual variance
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

  % seems that except for some prime m (e.g. 19), fft is always faster
  %p=500; m=20; a=rand(p,p,m);
  %tic; u=fft(a, m, 3); toc;
  %tic; F=fft(eye(m)).'; v=reshape(reshape(a,p*p,m)*F,p,p,m); toc;
  %max(abs((v-u)(:)))
end

% Copy from local repo test_blocklevinson
% Solve block toeplitz matrix by Levinson algorithm
% Derive from BlockLevinson_v3.m
% This version return mid results for constructing inverse of toeplitz matrix
% Ref: Xiao-Guang Lv and Ting-Zhu Huang(2013) The Inverses of Block Toeplitz Matrices.

function [fk, bk, ufk, ubk, A] = BlockLevinson_v3_inv(R, b_get_arcoef)
  if ~exist('b_get_arcoef','var')
    b_get_arcoef = false;
  end
  [p, m] = size(R);
  m = round(m/p)-1;

  % the 2d form [R(-m),...,R(-1),R(0)]. Note R(-k)=R(k)'
  b_R = reshape(flipdim(...
          permute(reshape(R,p,p,[]),[2,1,3])...
          , 3), p, []);

  y = [b_R(:,p+1:end-p) - R(:,p+1:end-p), zeros(p)];  % RHS for u tilde
  y1 = reshape(permute(flipdim(reshape(y,p,p,[]),3),[2,1,3]),p,[]).'; % RHS for u, (reverse order of y)
  if b_get_arcoef
    y = [y1, y', -R(:,p+1:end)'];  % add RHS for AR coef
  else
    y = [y1, y'];  % the RHS for u and u tilde transpose
  end
  T1 = R(:, 1:p);
  fk = inv(T1);
  bk = fk;
  xk = T1 \ y(1:p, :);
  p2 = size(y,2);

  for k=2:m
    eps_fk = b_R(:, end-k*p+1:end-p) * fk;
    eps_bk = R(:, p+1:k*p) * bk;
    Vfb = [[fk; zeros(p)], [zeros(p); bk]] / [eye(p), eps_bk; eps_fk, eye(p)];
    fk = Vfb(:,1:p);
    bk = Vfb(:,p+1:end);
    eps_xk = b_R(:, end-k*p+1:end-p) * xk(1:p*(k-1), :);
    xk = [xk;zeros(p,p2)] + bk * (y(p*(k-1)+1:p*k, :) - eps_xk);
  end

  % Break the 2d form (block column) to 3d form (block array)
  fk = permute(reshape(fk,p,m,[]),[1,3,2]);
  bk = permute(reshape(bk,p,m,[]),[3,1,2]); % add a transpose
  ufk = permute(reshape(xk(1:end,1:p),p,m,[]),[1,3,2]);
  ubk = permute(reshape(xk(1:end,p+1:p+p),p,m,[]),[3,1,2]); % add a transpose

  if b_get_arcoef
    A = xk(1:end,2*p+1:3*p)';
  end

end
