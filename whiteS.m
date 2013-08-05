% Help calculate GC from spectrum
% S is p*p*fftlen matrix

function S3 = whiteS(S)
if size(S,1)~=size(S,2)
  error('S shoule be p*p*fftlen matrix');
end

fftlen = size(S,3);
fq = (0:fftlen-1)/fftlen;

S11 = zeros(1,1,fftlen);
S11(1,1,:) = S(1,1,:);
%[Sn11,H11,de11] = sfactorization_wilson1(S11, fq);
[H11 de11] = S2H1D(S11);
S22 = zeros(1,1,fftlen);
S22(1,1,:) = S(2,2,:);
%[Sn22,H22,de22] = sfactorization_wilson1(S22, fq);
[H22 de22] = S2H1D(S22);

S3 = zeros(size(S));
for k = 1:fftlen
    iH = [inv(H11(1,1,k)),0;0, inv(H22(1,1,k))];
    S3(:,:,k) = iH*S(:,:,k)*iH';
end

