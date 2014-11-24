% Calculate spectrum according to AR coefficient
% A  : p * (p*m)      matrix, AR coefficient in 2-d form(no zeroth term)
% De : p * p          matrix, the residual covariance
% S  : p * p * fftlen  array, corresponding spectrum
% fqs: 1 * fftlen row vector, frequencies corresponding to the spectrum

function [S, fqs] = A2S(A, De, fftlen)

[p, m] = size(A);
m = round(m/p);
A = cat(3,eye(p),reshape(A,p,p,[]));  % convert A to 3-dim array
S = fft(A,fftlen,3);
hDe = chol(De,'lower');
for k = 1:fftlen
  h = S(:,:,k) \ hDe;
  S(:,:,k) = h*h';
end

fqs = ifftshift((0:fftlen-1)-floor(fftlen/2))'/fftlen;

end
