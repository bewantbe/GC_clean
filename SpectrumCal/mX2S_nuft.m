% Calculate spectrum from multi-trial non-uniformly sampled data mX
% mX is  p * len * n_trials  array, the sampled data
% mT is  len * n_trials  array, the sampling time
% aveS is fftlen*p*p matrix
% currently no window function applied
% does not subtract mean value from original data

function aveS = mX2S_nuft(mX, mT, fftlen)
mX = permute(mX,[2 1 3]);          % convert to len * p * n_trials
[len, p, n_trials] = size(mX);
if exist('fftlen','var') == 0
    fftlen = len;
end

% window function
wnd = ones(len,1);

desired_accuracy = 9;
aveS = zeros(fftlen,p,p);
S    = zeros(fftlen,p,p);
Jk   = zeros(fftlen, p);
% average over trials
for i_trial=1:n_trials
  % windowed Fourier transform
  for channel=1:p
    Jk(:,channel) = ifftshift(FGG_1d_type1(wnd .* mX(:,channel,i_trial), mT(:,i_trial), fftlen, desired_accuracy));
  end
  % get cross spectrum of one trial
  % due to symmetric of real data fft, this can be faster
  for chan1=1:p
    for chan2=1:p
      S(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
    end
  end
  aveS = aveS + S;
end
aveS = aveS / n_trials / fftlen;

