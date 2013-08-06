% Help calculate GC from spectrum
% S is p*p*fftlen matrix

function S3 = whiteS(S)
if size(S,1)~=size(S,2)
  error('S shoule be p*p*fftlen matrix');
end

[H11 de11] = S2H1D(squeeze(S(1,1,:)));
[H22 de22] = S2H1D(squeeze(S(2,2,:)));

S3 = zeros(size(S));
for k = 1:size(S,3)
    iH = [1/H11(k), 0; 0, 1/H22(k)];
    S3(:,:,k) = iH*S(:,:,k)*iH';
end
for k = 1:size(S,1)
    S3(k,k,:) = real(S3(k,k,:));
end
