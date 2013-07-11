% calculate GC approximately

function gc_app = getGCSapp(S_mt)

S3 = whiteS(S_mt);
de11 = mean(real(S3(1,1,:)));
de22 = mean(real(S3(2,2,:)));

covxy_f_app0 = real(ifft(squeeze(S3(1,2,:))'));

gc_app = zeros(2,2);
gc_app(2,1) = sum(covxy_f_app0(2:31).^2)/(de22*de11);
gc_app(1,2) = sum(covxy_f_app0(end-30:end).^2)/(de22*de11);

