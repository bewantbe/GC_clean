
function [aveS, fqs] = simple_spectrum(X, fftlen)
if ~exist('fftlen', 'var')
  fftlen = 1024;
  if size(X,2) < 3*fftlen
    fftlen = size(X,2)/3;
  end
end
fftslen = [fftlen, 0.5];       % 50% overlap
%f_wnd = @(x) 1 - 2*abs(x);   % bartlet window
f_wnd = @(x) 0.5+0.5*cos(2*pi*x);
[aveS, fqs] = mX2S_wnd(X, fftslen, f_wnd);
aveS = aveS(1:end/2,:,:);
fqs  = fqs(1:end/2);
