% Calculate spectrum by multitaper with DPSS (also called Slepian sequences)
% , which has max energy concentration in the main lobe
%
% aveS = mX2S_mt(mX, fftlen, nHBW [, n_taper])
%
% mX      : p * len * n_trials  array, the sampled data.
% fftlen  : Length of each small slice (like Bartlett's method).
%           pass [] means fftlen = len
% nHBW    : Half Band Width in unit of points in frequency domain.
% n_taper : Number of taper function, default is floor(nHBW/2-1).
% aveS    : len * p * p matrix, the averaged spectrum.
%
% Note: mean value of mX is NOT subtracted from the original data

function [aveS, fqs, winfo] = mX2S_mt(mX, fftlen, nHBW, n_taper)
[p, len, n_trials] = size(mX);  % conventional dimension order of xyy's code

if ~exist('fftlen','var') || isempty(fftlen)
  fftlen = len;
end
if length(fftlen)==1
  overlap_coef = 0;
else
  overlap_coef = fftlen(2);
  fftlen = fftlen(1);
  if overlap_coef>=1
    error('overlap_coef should smaller than 1');
  end
end

if ~exist('nHBW', 'var')
  nHBW = 3;
end
if ~exist('n_taper', 'var')
  n_taper = floor(2*nHBW-1);  % choose the tapers with fewer leakage
end
persistent tapers
persistent tapers_lambda
persistent weights
if size(tapers,1)~=fftlen || size(tapers,2)~=n_taper
  disp('re-generate dpss');
  [tapers, tapers_lambda] = dpss(fftlen, nHBW, n_taper);
  weights = tapers_lambda / sum(tapers_lambda);  % minor correction
  %weights = ones(n_taper, 1) / n_taper;
end

% break each trial into n_slice slices
n_slice = round((len-fftlen) / (fftlen*(1-overlap_coef))) + 1;
step_size = (len-fftlen) / max([n_slice-1, 1]); % not an integer in general

aveS = zeros(fftlen,p,p);
St   = zeros(fftlen,p,p);
% average over trials
for i_trial=1:n_trials
  for id_slice = 0:n_slice-1
    bg = round(id_slice*step_size);
    S = zeros(fftlen,p,p);
    % average over tapers
    for i_taper=1:n_taper
      % windowed Fourier transform
      Jk = fft( bsxfun(@times,...
        tapers(:,i_taper), mX(:,bg+1:bg+fftlen,i_trial).'...
        ));
      % Get cross spectrum of one slice
      for chan1=1:p
        St(:, chan1, chan1) = Jk(:,chan1).*conj(Jk(:,chan1));
        for chan2=chan1+1:p
          St(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
          St(:, chan2, chan1) = conj(St(:, chan1, chan2));
        end
      end
      S = S + weights(i_taper) * St;
    end
    aveS = aveS + S;
  end
end
aveS = aveS / n_trials / n_slice;
fqs = ifftshift((0:fftlen-1)-floor(fftlen/2))'/fftlen;

if nargout==3
  winfo.n_trials = n_trials;
  winfo.n_slices = n_slice * n_trials;
  winfo.n_overlap = fftlen - step_size;
  winfo.rel_error = 1/sqrt(n_trials*n_slice*n_taper);  % standard variance / amplitude
  fprintf('relative error(mt) = %.2f dB, band width = %f Hz\n',...
          10*log10(1+2*2*winfo.rel_error), 2*nHBW/fftlen);
end
