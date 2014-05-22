% calculate GC approximately
% S: spectrum fftlen*p*p
% ext_od: order used in the final summation

function [gc_app, de11, de22] = getGCSapp(S, ext_od)
if size(S,2)~=size(S,3)
  error('S shoule be fftlen*p*p matrix');
end
if ~exist('ext_od','var')
  ext_od = 50;
end
od = min(floor(size(S,1)/2), ext_od);

p = size(S,2);
WS = nStdWhiteS(S);
for k1=1:p
  gc_app(k1,k1) = 0;
  for k2=k1+1:p
    covxy_f = real(ifft(WS(:,k1,k2)));
    gc_app(k1,k2) = sum(covxy_f(2:od+1).^2);
    gc_app(k2,k1) = sum(covxy_f(end-od:end).^2);
  end
end
