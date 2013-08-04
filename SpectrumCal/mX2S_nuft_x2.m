% Calculate spectrum from multi-trial non-uniformly sampled data mX1, mX2
% mX1,mX2 is  p * len * n_trials  array, the sampled data
% mT1,mT2 is  len * n_trials  array, the sampling time
% aveS is fftlen*p*p matrix
% currently no window function applied
% does not subtract mean value from original data

function aveS = mX2S_nuft_x2(mX1, mX2, mT1, mT2, fftlen)
mX1 = permute(mX1,[2 1 3]);        % convert to len * p * n_trials
mX2 = permute(mX2,[2 1 3]);        % convert to len * p * n_trials
[len, p, n_trials] = size(mX1);
if exist('fftlen','var') == 0
    fftlen = len;
end

% window function
wnd = ones(len,1);

desired_accuracy = 6;
%nufft = @(X,T) ifftshift(FGG_1d_type1(fftshift(wnd .* X), fftshift(T), fftlen, desired_accuracy));

aveS = zeros(fftlen,p,p);
S    = zeros(fftlen,p,p);
Jk1  = zeros(fftlen, p);
Jk2  = zeros(fftlen, p);
% average over trials
for i_trial=1:n_trials
  % windowed Fourier transform
  for channel=1:p
    Jk1(:,channel) = ifftshift(FGG_1d_type1(wnd .* mX1(:,channel,i_trial), mT1(:,i_trial), fftlen, desired_accuracy));
    Jk2(:,channel) = ifftshift(FGG_1d_type1(wnd .* mX2(:,channel,i_trial), mT2(:,i_trial), fftlen, desired_accuracy));
    %Jk1(:,channel) = nufft(mX1(:,channel,i_trial), mT1(:,i_trial));
    %Jk2(:,channel) = nufft(mX2(:,channel,i_trial), mT2(:,i_trial));
  end
  % get cross spectrum of one trial
  % due to symmetric of real data fft, this can be faster
  for chan1=1:p
    for chan2=1:p
      S(:, chan1, chan2) = Jk1(:,chan1).*conj(Jk2(:,chan2));
    end
  end
  aveS = aveS + S;
end
aveS = aveS / n_trials / fftlen;

