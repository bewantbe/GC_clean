% calculate Self-Spectrum from X
% X is p*len matrix
% z is p*fftlen matrix
% no window function applied

function z = X2Sxx(X, fftlen)
[p, len] = size(X);
if exist('fftlen','var') == 0
    fftlen = 1024;
end
X = bsxfun(@minus, X, mean(X,2));  % must to strip means first
z = zeros(p, fftlen);
for ii=1:fftlen:len-fftlen+1
    y = fft(X(:, ii:ii+fftlen-1), fftlen, 2);
    z = z + y.*conj(y);
end
z = z / floor(len/fftlen) / fftlen;
end
