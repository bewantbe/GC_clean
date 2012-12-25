% Calculate covariance according to spectrum
% spectrum S is p * p * fftlen dim matrix
% m is the max offset of covariance R you want
% note that the length of fft should at least 2*(m+1), bigger can reduce bias problem

function R = S2cov(S, m)

[p, q, fftlen] = size(S);
covs = real(ifft(S,fftlen,3));
R = zeros(p,p*(m+1));
for k = 0 : m
	R(:,1+p*k:p+p*k) = covs(:,:,k+1);
end

end
