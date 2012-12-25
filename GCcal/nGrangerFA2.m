% Granger Causality in frequency domain. only a prototype

function [gcy2x, gcx2y, fx2y, fy2x, Sxx, Syy, Hw] = nGrangerFA2(A, noisecov, fftlen)
if (exist('fftlen', 'var')==0)
    fftlen = 1024;
end

P = [1 0; -noisecov(1,2)/noisecov(1,1) 1];
Q = [1 -noisecov(1,2)/noisecov(2,2); 0 1];

[p, m] = size(A);
m = round(m/p);

vA = zeros(p,p,m+1);
vA(:,:,1) = eye(p);             % convert A to 3-dim array
for k = 1 : m
    vA(:,:,k+1) = A(:,1-p+p*k:p*k);
end
vAw = fft(vA,fftlen,3);
Hw  = zeros(p,p,m);              % transfer function
pHw = zeros(p,p,m);              % transfer function
qHw = zeros(p,p,m);              % transfer function
S   = zeros(p,p,m);
pS  = zeros(p,p,m);
qS  = zeros(p,p,m);
pnoisecov = P * noisecov * P';
qnoisecov = Q * noisecov * Q';
for k = 1 : fftlen
    Hw(:,:,k) = inv(vAw(:,:,k));
    pHw(:,:,k) = inv(P * vAw(:,:,k));
    qHw(:,:,k) = inv(Q * vAw(:,:,k));
    S (:,:,k) = Hw(:,:,k) * noisecov * Hw(:,:,k)';
    pS(:,:,k) = pHw(:,:,k) * pnoisecov * pHw(:,:,k)';
    qS(:,:,k) = qHw(:,:,k) * qnoisecov * qHw(:,:,k)';
end

Sxx = real(squeeze(S(1,1,:)));
Syy = real(squeeze(S(2,2,:)));

s_fq = 1:fftlen/2;

pHxx = squeeze(pHw(1,1,:));
qHyy = squeeze(qHw(2,2,:));

pSx0 = abs(pHxx.*conj(pHxx))*noisecov(1,1);
qSy0 = abs(qHyy.*conj(qHyy))*noisecov(2,2);

fy2x = log(Sxx./pSx0);
fx2y = log(Syy./qSy0);
%{
figure(1);
set(gca, 'fontsize', 15);
set(gca, 'fontweight', 'bold');
plot(s_fq/fftlen, [Sxx(s_fq)'; Syy(s_fq)'], 'linewidth', 2);
legend('Sxx','Syy');
print('-dpng','pic/s.png','-r100');

figure(2);
set(gca, 'fontsize', 15);
set(gca, 'fontweight', 'bold');
plot(s_fq/fftlen, [fy2x(s_fq)'; fx2y(s_fq)']*1000, 'linewidth', 2);
%ax=axis();
%ax(4)=1.2;
%axis(ax);
legend('f(y -> x) * 1000','f(x -> y) * 1000');
print('-dpng','pic/fw.png','-r100');

figure(3);
set(gca, 'fontsize', 15);
set(gca, 'fontweight', 'bold');
plot(s_fq/fftlen, qSy0(s_fq), 'linewidth', 2);
legend('qHyy*noisecov(2,2)*qHyy^T');
print('-dpng','pic/sy.png','-r100');

figure(4);
set(gca, 'fontsize', 15);
set(gca, 'fontweight', 'bold');
plot(s_fq/fftlen, (Syy(s_fq) - qSy0(s_fq))*1000, 'linewidth', 2);
legend('(Syy(s_fq) - qHyy*noisecov(2,2)*qHyy^T) * 1000');
print('-dpng','pic/sy-sy0.png','-r100');
%}
gcy2x = mean(fy2x);
gcx2y = mean(fx2y);

end
