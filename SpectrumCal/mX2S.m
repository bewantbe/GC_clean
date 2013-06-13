% Calculate Spectrum from multi-trial non-uniformly sampled data mX
% mX is  p * len * n_trials  array, the sampled data
% (optional) mT is  len * n_trials  array, the sampling time
% currently no window function applied
% does not subtract mean value from original data

function aveS = mX2S(mX, mT, fftlen)
mX = permute(mX,[2 1 3]);
[len, p, n_trials] = size(mX);
if exist('mT','var') && length(mT)==1   % mT is scaler, called as if mX2S(mX, fftlen)
  fftlen = mT;
  clear('mT');
end
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

%if exist('mT','var')
%  % do NUFT
%  desired_accuracy = 9;
%  ft = @(X,T) ifftshift(FGG_1d_type1(fftshift(wnd .* X), fftshift(T), fftlen, desired_accuracy));
%  stft = 'ft(mX(:,channel,i_trial), mT(:,i_trial))';
%else
%  ft = @(X) fft(wnd .* X, fftlen);
%  stft = 'ft(mX(:,channel,i_trial))';
%end

if exist('mT','var')
  % do NUFT
  desired_accuracy = 9;
  ft = @(channel,i_trial) ifftshift(FGG_1d_type1(fftshift(wnd .* mX(:,channel,i_trial)), fftshift(mT(:,i_trial)), fftlen, desired_accuracy));
else
  ft = @(X) fft(wnd .* X, fftlen);
  stft = 'ft(mX(:,channel,i_trial))';
end

aveS = zeros(fftlen, p,p);
S    = zeros(fftlen, p,p);
% average over trials
for i_trial=1:n_trials
  % windowed Fourier transform
  Jk = zeros(fftlen, p);
  for channel=1:p
    Jk(:,channel) = eval(stft);
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

