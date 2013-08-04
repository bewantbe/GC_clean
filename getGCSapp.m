% calculate GC approximately
% S_mt: fftlen*p*p

function gc_app = getGCSapp(S_mt)
S_mt = ipermute(S_mt, [3,1,2]);   % convert to p*p*fftlen
if size(S_mt,1)~=size(S_mt,2)
  error('S shoule be p*p*fftlen matrix');
end

S3 = whiteS(S_mt);
de11 = mean(real(S3(1,1,:)));
de22 = mean(real(S3(2,2,:)));

covxy_f_app0 = real(ifft(squeeze(S3(1,2,:))'));

gc_app = zeros(2,2);
gc_app(2,1) = sum(covxy_f_app0(2:31).^2)/(de22*de11);
gc_app(1,2) = sum(covxy_f_app0(end-30:end).^2)/(de22*de11);

