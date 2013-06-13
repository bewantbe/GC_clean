% Calculate spectrum by multitaper with DPSS
% Result is the averaged spectrum by trials in the given data
% input: size(data) = #variables * #trials * #datapoints_in_each_experiment
% does not subtract mean value from original data

function aveS = mX2S_mt(mX, fftlen, halfbandwidth, n_taper)

% for compatibility of origninal non-parametric code
%[p, n_trials, len] = size(mX);
%mX = permute(mX,[3 1 2]);       % it's convenient to use this dimension order below
%[len, p, n_trials] = size(mX);

[p, len, n_trials] = size(mX);  % conventional dimension order of xyy's code
mX = permute(mX,[2 1 3]);       % it's convenient to use this dimension order below
%[len, p, n_trials] = size(mX);

if ~exist('fftlen','var')
  % which one is proper?
%  fftlen = 2^(nextpow2(len)+1);
  fftlen = len;
end

% Use DPSS tapers (also called Slepian sequences)
% in order to maximize the energy concentration in the main lobe
if ~exist('halfbandwidth', 'var')
  halfbandwidth = 3;
  %n_taper = 8;                       % number of tapers
end
if ~exist('n_taper', 'var')
  tapers = dpss(len, halfbandwidth);
else
  tapers = dpss(len, halfbandwidth, n_taper);
end

aveS = zeros(p,p,fftlen);
St   = zeros(p,p,fftlen);
% average over trials
for i_trial=1:n_trials
  % average over tapers
  S = zeros(p,p,fftlen);
  for i_taper=1:n_taper
    % windowed Fourier transform
    Jk = zeros(fftlen, p);
    for channel=1:p
      %Jk(:,channel) = fft(tapers(:,i_taper) .* (mX(:,channel,i_trial) - mean(mX(:,channel,i_trial))), fftlen);
      Jk(:,channel) = fft(tapers(:,i_taper) .* mX(:,channel,i_trial), fftlen);
    end
    % get cross spectrum of one tapered data
    % due to symmetric of real data fft, this can be faster
    for chan1=1:p
      for chan2=1:p
        St(chan1, chan2, :) = Jk(:,chan1).*conj(Jk(:,chan2));
      end
    end
    S = S + St;
  end
  aveS = aveS + S / n_taper;
end
aveS = aveS / n_trials;

