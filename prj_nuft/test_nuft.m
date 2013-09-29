% test spectrum

elen = 1000;
simu_time = 1e6;
De = [1.0, 0.4; 0.4, 0.7];
A = [-0.9 ,  0.0, 0.5, 0.0;
     -0.16, -0.8, 0.2, 0.5];
#S_ans = A2S(A, De, 1000);  S_ans = permute(S_ans,[3,1,2]);

tic();
oX = gendata_linear(A, De, simu_time+elen);
oX(:,1:elen) = [];
#oX = gdata(simu_time, 2, 2);
toc();

resample_mode = 'r';
T_segment = 1e3;
stv0 = 1;
stv  = 1;

tic();
mlen = T_segment/stv0;
slen = round(T_segment/stv);
[X1,T1] = SampleNonUnif(oX, mlen, slen, resample_mode,
                        simu_time/T_segment*stv/stv0);
#[X1,T1] = SampleNonUnif(oX, mlen, slen, resample_mode);
toc();

tic();
fftlen = mlen;
[S1, fqs] = mX2S_nuft(X1, (T1-1)/mlen, fftlen,
                      simu_time/T_segment*stv/stv0);
fqs = fqs/T_segment;
toc();

tic();
[X0,T0] = SampleNonUnif(oX, mlen, slen, 'u');
[Su, fqsu] = mX2S_ft(X0);
fqsu = fqsu/T_segment;
toc();

DeX = getcovpd(oX, 0);
#S1 = nuft_bias_removal(S1, DeX, stv, slen);
S1(:,1,1) = S1(:,1,1) - DeX(1,1);
S1(:,2,2) = S1(:,2,2) - DeX(2,2);

figure(1);
plot(fftshift(fqs), fftshift(S1(:,1,1)), fftshift(fqsu), fftshift(Su(:,1,1)));
#plot(fftshift(fqs), fftshift(S1(:,1,1)) - fftshift(Su(:,1,1)));
#mean(real(fftshift(S1(:,1,1))-fftshift(Su(:,1,1))))
#mean(real(fftshift(S1(:,2,2))-fftshift(Su(:,2,2))))

