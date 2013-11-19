%Fast nufft, up to frequency fq_max
% T should in range [0, 2*pi]
% usually fqs(k) = -M/2 .. M/2-1
% X is p*len
% borrow from test_nufft/nufftw.m

function [F_r, fqs] = nufftw(X, T, fq_max)
if any(T(:)>2*pi) || any(T(:)<0)
  error('T must in the range [0, 2*pi]');
end

along_dim = 1;
if (size(X,1)<size(X,2))
  X = X.';
  T = T.';
  along_dim = 2;
end

N = size(X,1);                  % length of data points
M = round(2*fq_max);            % number of frequency points
fqs = (1:M) - floor(M/2) - 1;   % desired frequency points
fqs = ifftshift(fqs)';          % so consist with fft()

R = 2;
M_r = M*R;
%M_r = lowest_smooth_number_fast(M_r);  % choose a number that only have
M_r = lowest_smooth_number_exact(M_r);  % choose a number that only have
R   = M_r/M;                           %  factor 2, 3, 5. so fft faster
%fprintf('use M_r = %d\n', M_r);
M_sp= 16;

tau = pi * M_sp / M^2 / (R*(R-0.5));

f_r = GaussianConvGrid(X, T, M_r, M_sp, tau);

F_r = fft(f_r)/M_r;

F_r = F_r([1:ceil(fq_max-0.25), M_r-floor(fq_max-0.75):M_r], :);

G = sqrt(pi/tau) * exp(fqs.^2 * tau);
F_r = bsxfun(@times, F_r, G);

if (along_dim == 2)
  F_r = F_r.';
  fqs = fqs';
end
