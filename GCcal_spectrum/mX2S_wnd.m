% Calculate spectrum from multi-trial uniformly sampled data mX
%
% [aveS, fqs] = mX2S_wnd(mX, fftlen, f_wnd)
%
% mX     : p * len * n_trials  array, the sampled data
% fftlen : Length of each small slice (Bartlett's method).
%          Possible to specify the portion of overlap as second component (Welch's method).
% f_wnd  : The window function, either a vector that have the same len as mX
%          or a function defined in (-0.5,0.5).
% aveS   : fftlen * p * p  matrix, the averaged spectrum
%
% Note: mean value of mX is NOT subtracted from the original data.
%
% Usage Example:
%
%  eX = sin(cumsum(0.4+0.05*randn(1, 1e6)));
%  fftlen = [1000, 0.5];       % 50% overlap
%  f_wnd = @(x) 1 - 2*abs(x);  % bartlet window
%  [aveS, fqs] = mX2S_wnd(eX, fftlen, f_wnd);
%  aveS = fftshift(aveS);      % convenient for plot
%  fqs  = fftshift(fqs);
%  semilogy(fqs, aveS);

function [aveS, fqs, winfo] = mX2S_wnd(mX, fftlen, f_wnd)
[p, len, n_trials] = size(mX);
if exist('fftlen','var') == 0 || isempty(fftlen)
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
fqs = ifftshift((0:fftlen-1)-floor(fftlen/2))'/fftlen;

% get window function
wnd = ones(fftlen, 1);      % use column vector
if exist('f_wnd','var') && ~isempty(f_wnd)
  if ismatrix(f_wnd) && isnumeric(f_wnd)
    wnd(:) = f_wnd;
  else
    % f_wnd is defined in [-0.5, 0.5]
    wnd(:) = f_wnd( ((1:fftlen)-0.5)/fftlen - 0.5 );
  end
end

% break each trial into n_slice slices
n_slice = round((len-fftlen) / (fftlen*(1-overlap_coef))) + 1;
step_size = (len-fftlen) / max([n_slice-1, 1]); % not an integer in general

aveS = zeros(fftlen,p,p);
S    = zeros(fftlen,p,p);
% average over trials
for i_trial=1:n_trials
  for id_slice = 0:n_slice-1
    bg = round(id_slice*step_size);
    % windowed Fourier transform
    Jk = fft(bsxfun(@times, wnd, mX(:,bg+1:bg+fftlen,i_trial).'));
    % get cross spectrum of one slice
    for chan1=1:p
      S(:, chan1, chan1) = Jk(:,chan1).*conj(Jk(:,chan1));
      for chan2=chan1+1:p
        S(:, chan1, chan2) = Jk(:,chan1).*conj(Jk(:,chan2));
        S(:, chan2, chan1) = conj(S(:, chan1, chan2));
      end
    end
    aveS = aveS + S;
  end
end
aveS = aveS / (n_trials * n_slice * (wnd'*wnd));

if nargout==3
  winfo.n_trials = n_trials;
  winfo.n_slices = n_slice * n_trials;
  winfo.n_overlap = fftlen - step_size;
  winfo.rel_error = 1/sqrt(n_trials*n_slice); % standard variance / amplitude
  fprintf('relative error = %.2f dB\n', 10*log10(1+2*2*winfo.rel_error));
end
