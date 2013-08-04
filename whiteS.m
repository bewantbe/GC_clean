% Help calculate GC from spectrum
% S_mt is p*p*fftlen matrix

function S3 = whiteS(S_mt)
if size(S_mt,1)~=size(S_mt,2)
  error('S shoule be p*p*fftlen matrix');
end

fftlen = size(S_mt,3);
fq = (0:fftlen-1)/fftlen;

S11 = zeros(1,1,fftlen);
S11(1,1,:) = S_mt(1,1,:);
[Sn11,H11,de11] = sfactorization_wilson1(S11, fq);
S22 = zeros(1,1,fftlen);
S22(1,1,:) = S_mt(2,2,:);
[Sn22,H22,de22] = sfactorization_wilson1(S22, fq);

S3 = zeros(size(S_mt));
for k = 1:fftlen
    iH = [inv(H11(1,1,k)),0;0, inv(H22(1,1,k))];
    S3(:,:,k) = iH*S_mt(:,:,k)*iH';
end

