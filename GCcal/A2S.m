% Calculate spectrum according to AR coefficient A
% A is p * (p*m) dim matrix
% noisecov is p * p matrix, represents the variance of noise term
% S is p * p * fftlen dim matrix

function S = A2S(A, noisecov, fftlen)

[p, m] = size(A);
m = round(m/p);
vA(:,:,1) = eye(p);             % convert A to 3-dim array
for k = 1 : m
    vA(:,:,k+1) = A(:,1-p+p*k:p*k);
end
vAw = fft(vA,fftlen,3);
Hw = zeros(p,p,m);              % transfer function
for k = 1 : fftlen
    Hw(:,:,k) = inv(vAw(:,:,k));
    vAw(:,:,k) = Hw(:,:,k) * noisecov * Hw(:,:,k)';
end
S = vAw;

end
