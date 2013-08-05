% calculate GC approximately
% S: fftlen*p*p

function [gc_app, de11, de22] = getGCSapp(S)
S = ipermute(S, [3,1,2]);   % convert to p*p*fftlen
if size(S,1)~=size(S,2)
  error('S shoule be p*p*fftlen matrix');
end

S3 = whiteS(S);
de11 = mean(real(S3(1,1,:)));
de22 = mean(real(S3(2,2,:)));

covxy_f_app0 = real(ifft(squeeze(S3(1,2,:))'));

od = min(floor(size(S,1)/2),30);
gc_app = zeros(2,2);
gc_app(2,1) = sum(covxy_f_app0(2:od+1).^2)/(de22*de11);
gc_app(1,2) = sum(covxy_f_app0(end-od:end).^2)/(de22*de11);

