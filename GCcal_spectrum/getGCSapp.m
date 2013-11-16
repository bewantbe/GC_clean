% calculate GC approximately
% S: spectrum fftlen*p*p
% ext_od: order used in the final summation

function [gc_app, de11, de22] = getGCSapp(S, ext_od)
if size(S,2)~=size(S,3)
  error('S shoule be fftlen*p*p matrix');
end
if ~exist('ext_od','var')
  ext_od = 30;
end

[S3, de11, de22] = StdWhiteS(S);
covxy_f_app0 = real(ifft(S3(:,1,2)));
od = min(floor(size(S3,1)/2), ext_od);
gc_app = zeros(2,2);
gc_app(1,2) = sum(covxy_f_app0(2:od+1).^2);
gc_app(2,1) = sum(covxy_f_app0(end-od:end).^2);
