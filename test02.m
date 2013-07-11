%

X = gendata_neu('net_2_2', 0.02, 1, 0.012, 1e6, 0.5);
X = bsxfun(@minus, X, mean(X,2));

mlen = 512;                                % length of each trials at most
n_trials = floor(size(X,2)/mlen);          % number of trials at most
mX = reshape(X(:,1:mlen*n_trials), size(X,1), mlen, []);

S_mt = mX2S_mt(mX);

fftlen = size(S_mt,3);
fq = (0:fftlen-1)/fftlen;
[Snew,H,de] = sfactorization_wilson1(S_mt, fq);

cfq = fq - 0.5;

%figure(1);
%semilogy(cfq, fftshift(squeeze(S_mt(1,1,:))),  cfq, fftshift(squeeze(Snew(1,1,:))));
%xlim([-0.5,0.5]);


S2 = zeros(size(S_mt));
for k = 1:fftlen
%    iH = inv(H(:,:,k));
    iH = inv([H(1,1,k),0;0,H(2,2,k)]);
    S2(:,:,k) = iH*S_mt(:,:,k)*iH';
end


S11 = zeros(1,1,fftlen);
S11(1,1,:) = S_mt(1,1,:);
[Sn11,H11,de11] = sfactorization_wilson1(S11, fq);
S22 = zeros(1,1,fftlen);
S22(1,1,:) = S_mt(2,2,:);
[Sn22,H22,de22] = sfactorization_wilson1(S22, fq);

S3 = zeros(size(S_mt));

for k = 1:fftlen
    iH = [inv(H11(1,1,k)),0;0, inv(H22(1,1,k))];
    S3(:,:,k) = iH*S_mt(:,:,k)*iH';
end


figure(6);
plot(
cfq, fftshift(real(squeeze(S2(1,1,:)))),
cfq, fftshift(real(squeeze(S2(2,2,:)))),
cfq, fftshift(real(squeeze(S3(1,1,:)))),
cfq, fftshift(real(squeeze(S3(2,2,:))))
);
xlim([-0.5,0.5]);

figure(7);
plot(
cfq, fftshift(real(squeeze(S2(1,2,:)))),
cfq, fftshift(imag(squeeze(S2(1,2,:)))),
cfq, fftshift(real(squeeze(S3(1,2,:)))),
cfq, fftshift(imag(squeeze(S3(1,2,:))))
);
legend('Re(S2)', 'Im(S2)', 'Re(S3)', 'Im(S3)');
xlim([-0.5,0.5]);

covxy_f_app0 = real(ifft(squeeze(S2(1,2,:))'));
figure(8);
covxy_f_app = fftshift(covxy_f_app0);
%rg = (abs(cfq)<0.1) & (cfq ~= 0);
rg = (abs(cfq)<0.1);
plot(cfq(rg), covxy_f_app(rg), '-o');

% approximate GC
disp('GC x->y');
sum(covxy_f_app0(1:20).^2)/(de22*de11)
disp('GC y->x');
sum(covxy_f_app0(end-20:end).^2)/(de22*de11)



 SGrangerS(S3)
