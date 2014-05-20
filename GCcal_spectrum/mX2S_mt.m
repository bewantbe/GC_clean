% Calculate spectrum by multitaper with DPSS (also called Slepian sequences)
% , which has max energy concentration in the main lobe
%
% aveS = mX2S_mt(mX, fftlen, nHBW [, n_taper])
%
% mX      : p * len * n_trials  array, the sampled data.
% fftlen  : Should be the same as len,you may pass an empty matrix [] to it.
% nHBW    : Half Band Width in unit of points in frequency domain.
% n_taper : Number of taper function, default is floor(nHBW/2).
% aveS    : len * p * p matrix, the averaged spectrum.
%
% Note: mean value of mX is NOT subtracted from the original data

function [aveS, fqs, winfo] = mX2S_mt(mX, fftlen, nHBW, n_taper)
[p, len, n_trials] = size(mX);  % conventional dimension order of xyy's code
mX = permute(mX,[2 1 3]);       % convert to len * p * n_trials, for speed

if ~exist('fftlen','var') || isempty(fftlen)
  fftlen = len;
end

if ~exist('nHBW', 'var')
  nHBW = 3;
end
if ~exist('n_taper', 'var')
  n_taper = floor(2*nHBW-1);  % choose the tapers with fewer leakage
end
persistent tapers
if size(tapers,1)~=fftlen || size(tapers,2)~=n_taper
  disp('re-generate dpss');
  tapers = dpss(fftlen, nHBW, n_taper);
end

aveS = zeros(fftlen,p,p);
St   = zeros(fftlen,p,p);
Jk   = zeros(fftlen, p);
% average over trials
for i_trial=1:n_trials
  S = zeros(fftlen,p,p);
  % average over tapers
  for i_taper=1:n_taper
    % windowed Fourier transform
    for channel=1:p
      Jk(:,channel) = fft(tapers(:,i_taper) .* mX(:,channel,i_trial), fftlen);
    end
    % get cross spectrum of one tapered data
    % due to symmetric of real data fft, this might be faster
    for chan1=1:p
      for chan2=1:p
        St(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
      end
    end
    S = S + St;
  end
  aveS = aveS + S / n_taper;
end
aveS = aveS / n_trials;
fqs = ifftshift((0:fftlen-1)-floor(fftlen/2))'/fftlen;

if nargout==3
  winfo.rel_error = 1/sqrt(n_trials*n_taper);  % standard variance / amplitude
  fprintf('relative error(mt) = %.2f dB\n', 10*log10(1+2*2*winfo.rel_error));
end
