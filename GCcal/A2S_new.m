% Calculate spectrum according to AR coefficient A
% A is p * (p*m) dim matrix
% De is p * p matrix, represents the variance of noise term
% S is p * p * fftlen dim matrix

function [S fqs] = A2S_new(A, De, fftlen)

[p, m] = size(A);
m = round(m/p);
A = cat(3,eye(p),reshape(A,p,p,[]));  % convert A to 3-dim array
G = fft(A,fftlen,3);
S = zeros(p,p,m);
hDe = chol(De,'lower');
for k = 1:fftlen
  h = G(:,:,k) \ hDe;
  S(:,:,k) = h*h';
end

fqs = ifftshift((0:fftlen-1)-floor(fftlen/2))'/fftlen;

end
