% Calculate Spectrum from multi-trial uniformly sampled data mX
% mX is  len * p * n_trials  array
% currently no window function applied
% does not subtract mean value from original data

function aveS = mX2S_ft(mX, fftlen)
mX = permute(mX,[2 1 3]);
[len, p, n_trials] = size(mX);
if exist('fftlen','var') == 0
    fftlen = len;
end

% window function
wnd = ones(len,1);
% Blackman window, am I correctly use it?
%w0  = (2*pi*(0:len-1)/(len-1))';
%wnd = 7938/18608 - 9240/18608*cos(w0) + 1430/18608*cos(2*w0);
% Hanning window
%w0  = (2*pi*(0:len-1)/(len-1))';
%wnd = sin(w0).^2;

aveS = zeros(p,p,fftlen);
S    = zeros(p,p,fftlen);
% average over trials
for i_trial=1:n_trials
  % windowed Fourier transform
  Jk = zeros(fftlen, p);
  for channel=1:p
%    Jk(:,channel) = fft(wnd .* (mX(:,channel,i_trial) - mean(mX(:,channel,i_trial))), fftlen);
    Jk(:,channel) = fft(wnd .* mX(:,channel,i_trial), fftlen);
  end
  % get cross spectrum of one trial
  % due to symmetric of real data fft, this can be faster
  for chan1=1:p
    for chan2=1:p
      S(chan1, chan2, :) = Jk(:,chan1).*conj(Jk(:,chan2));
    end
  end
  aveS = aveS + S;
end
aveS = aveS / n_trials / fftlen;

