% Calculate spectrum from multi-trial non-uniformly sampled data mX
% Input:
%   mX is  p * len * n_trials  array, the sampled data
%   mT is  len * n_trials  array, the sampling time, in range [0,1]
%      or  p * len * n_trials  array, means all time point can be different
% Output:
%   aveS is fftlen*p*p matrix
%   correcponding frequencies: fqs = 0, 1,..., M/2-1, -M/2, -M/2+1,...,-1
% Currently no window function applied
% CAUTION: this function does not subtract mean value from original data

function [aveS, fqs] = mX2S_nuft(mX, mT, fftlen, TT)
mX = permute(mX,[2 1 3]);        % convert to len * p * n_trials
[len, p, n_trials] = size(mX);
b_uniform_time = true;
if n_trials == 1 && ndims(mT) == 2 || ndims(mT) == 3
    b_uniform_time = false;
    mT = permute(mT,[2 1 3]);
end
if exist('fftlen','var') == 0
    fftlen = len;
end
fqs = ifftshift((0:fftlen-1)-floor(fftlen/2));

aveS = zeros(fftlen,p,p);
S    = zeros(fftlen,p,p);
Jk   = zeros(fftlen, p);
% average over trials
for i_trial=1:n_trials
  if b_uniform_time
    for channel=1:p
      Jk(:,channel) = nufftw(mX(:,channel,i_trial),...
                             2*pi*mT(:,i_trial), fftlen/2);
    end
  else
    for channel=1:p
      Jk(:,channel) = nufftw(mX(:,channel,i_trial),...
                             2*pi*mT(:,channel,i_trial), fftlen/2);
    end
  end
  % get cross spectrum of one trial
  % due to symmetric of real data fft, this could be faster
  for chan1=1:p
    for chan2=1:p
      S(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
    end
  end
  aveS = aveS + S;
end
if ~exist('TT','var')
  % scale data according to input data length, instead of fftlen
  aveS = aveS / n_trials / len;
else
  % normalize to power spectral density of continuous time stationary process
  aveS = aveS*TT / n_trials / len^2;
end
