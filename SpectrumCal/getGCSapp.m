% calculate GC approximately
% S: fftlen*p*p

function [gc_app, de11, de22] = getGCSapp(S)
if size(S,2)~=size(S,3)
  error('S shoule be fftlen*p*p matrix');
end

%S = ipermute(S, [3,1,2]);   % convert to p*p*fftlen
%S3 = whiteS(S);
%de11 = mean(real(S3(1,1,:)));
%de22 = mean(real(S3(2,2,:)));

%covxy_f_app0 = real(ifft(squeeze(S3(1,2,:))'));

%od = min(floor(size(S,1)/2),30);
%gc_app = zeros(2,2);
%gc_app(2,1) = sum(covxy_f_app0(2:od+1).^2)/(de22*de11);
%gc_app(1,2) = sum(covxy_f_app0(end-od:end).^2)/(de22*de11);

[S3, de11, de22] = StdWhiteS(S);
covxy_f_app0 = real(ifft(S3(:,1,2)));
od = min(floor(size(S3,1)/2),30);
gc_app = zeros(2,2);
gc_app(2,1) = sum(covxy_f_app0(2:od+1).^2);
gc_app(1,2) = sum(covxy_f_app0(end-od:end).^2);
