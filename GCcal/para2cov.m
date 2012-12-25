% Calculate cov though AR parameters directly, v1.0
% Using spectrum to calculate order m, p-var AR model covariance series
% extm means extra order m to calculate (covariance series)

% only support following calls:
%  para2cov(1, 2)
%  para2cov(2, 2)
%  para2cov(2, 3)
%  para2cov(3, 5)

function R = para2cov(m, p, extm)
if (exist('m','var')==0 || exist('p','var')==0)
	error('usage: para2cov(2, 2) or para2cov(2, 3) or para2cov(3, 5)');
end
if (exist('extm','var')==0)
	extm = 0;
end
if (m==1 && p==2)
	noisecov = [1.0, 0.5; 0.5, 1.0];
	A = [ 0.6, -0.3;
	     -0.2,  0.8];
end
if (m==2 && p==2)
	noisecov = [1.0, 0.4; 0.4, 0.7];
	A = [-0.9 ,  0.0, 0.5, 0.0;
	     -0.16, -0.8, 0.2, 0.5];
end
if (m==2 && p==3)
	noisecov = diag([0.3 1.0 0.2]);
	A = [-0.8  0.0 -0.4  0.5 -0.2  0.0;
	      0.0 -0.9  0.0  0.0  0.8  0.0;
	      0.0 -0.5 -0.5  0.0  0.0  0.2];
end
if (m==3 && p==5)
	noisecov = diag([0.6 0.5 0.3 0.3 0.6]);
	s2 = sqrt(2);
	A = zeros(5, 5*3);
	A(1,1) =-0.95*s2;  A(1,6) = 0.9025;
	A(2,6) =-0.5;
	A(3,11)= 0.4;
	A(4,6) = 0.5;      A(4,4) =-0.25*s2;  A(4,5) =-0.25*s2;
	A(5,4) = 0.25*s2;  A(5,5) =-0.25*s2;
end

if (exist('A','var')==0 || exist('noisecov','var')==0)
	error('Not support these AR parameters.');
end

vA(:,:,1) = eye(p);             % convert A to 3-dim array
for k = 1 : m
	vA(:,:,k+1) = A(:,1-p+p*k:p*k);
end
fftlen = 1024;                  % at least 2*(m+1), bigger can reduce bias problem
vAw = fft(vA,fftlen,3);
Hw = zeros(p,p,m);              % transfer function
for k = 1 : fftlen
	Hw(:,:,k) = inv(vAw(:,:,k));
	vAw(:,:,k) = Hw(:,:,k) * noisecov * Hw(:,:,k)';
end
covs = real(ifft(vAw,fftlen,3));
m = m + extm;                     % want extra covariance?
R = zeros(p,p*(m+1));             % return only first m covariance matrices
for k = 0 : m
	R(:,1+p*k:p+p*k) = covs(:,:,k+1);
end
end
