% Help calculate GC from spectrum
% S is p*p*fftlen matrix

function S3 = StdWhiteS(S)
if size(S,2)~=size(S,3)
  error('S shoule be fftlen*p*p matrix');
end

X = S2X1D(S(:,1,1));
Y = S2X1D(S(:,2,2));

S3 = zeros(size(S));
%S3(:,1,1) = 1;
%S3(:,1,2) = (1./X) .* squeeze(S(:,1,2)) .* conj(1./Y);
%S3(:,2,1) = (1./Y) .* squeeze(S(:,2,1)) .* conj(1./X);
%S3(:,2,2) = 1;

S3(:,1,1) = (1./X) .* S(:,1,1) .* conj(1./X);
S3(:,1,2) = (1./X) .* S(:,1,2) .* conj(1./Y);
S3(:,2,1) = (1./Y) .* S(:,2,1) .* conj(1./X);
S3(:,2,2) = (1./Y) .* S(:,2,2) .* conj(1./Y);
