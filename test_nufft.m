% test nufft
tocs = @(st) fprintf('%s: t = %6.3fs\n', st, toc());

stv0 = 0.125;
oX = gendata_neu('net_2_2', 0.02, 1, 0.012, 1e5, stv0);
oX = bsxfun(@minus, oX, mean(oX,2));
p = size(oX,1);

T_segment = 0.5*128;                  % ms

slen = 128;                           % points in one non-unif sample

% Unif
stv = 4*stv0;
uX = oX(:,1:(stv/stv0):end);
mlen = T_segment/stv;                 % length of each trials
n_trials = floor(size(uX,2)/mlen);    % number of trials at most
mX = reshape(uX(:,1:mlen*n_trials), size(uX,1), mlen, []);
ufq = (0:mlen-1)/mlen/stv;

tic;
S_us = mX2S_ft(mX);
tocs('ft      ');
S_us = S_us*slen/mlen;

% non-Unif
rand('state',2);
tic;
mlen = T_segment/stv0;
[X1,X2,T1,T2] = SampleNonUnif(oX, mlen, slen, '1');
fq  = (0:slen-1)/T_segment;
fqs = [0:slen/2-1, -slen/2:-1]/T_segment;
fq_max = slen/T_segment;
tocs('resample');

tic;
S_nu = mX2S_nuft(X1, T1);
tocs('nufft   ');

tic;
S_x2 = mX2S_nuft_x2(X1, X2, T1, T2);
tocs('nufft_x2');
S_x2_bk = S_x2;

S1 = Makeup4SpectrumFact(S_x2_bk);

%S_x2(:,1,1) = (S_x2(:,1,1)+fliplr(S_x2([2:slen,1],1,1)))/2;
%S_x2(:,2,2) = (S_x2(:,2,2)+fliplr(S_x2([2:slen,1],2,2)))/2;
%S_x2(:,1,1) = abs(S_x2(:,1,1));
%S_x2(:,2,2) = abs(S_x2(:,2,2));

v = var(oX,[],2);
% plot
figure(1);
%plot(fq, real(S_nu(:,1,1))-v(1), fq, abs(S_x2(:,1,1)), ufq, real(S_us(:,1,1)), '-o');
%legend('nufft-bia','nufft-x2','unif');  xlim([0,fq_max]);
plot(fq, real(S1(:,1,1)), fq, real(S_x2(:,1,1)), ufq, real(S_us(:,1,1)), '-o');
legend('nufft-S1','nufft-x2','unif');  xlim([0,fq_max]);

%figure(2);
%plot(fq, angle(S_x2(:,1,1)), '-o');  xlim([0,fq_max]);

figure(3);
plot(fq, abs(S_nu(:,1,2)), fq, abs(S_x2(:,1,2)), ufq, abs(S_us(:,1,2)), '-o');
legend('nufft-bia','nufft-x2','unif');  xlim([0,fq_max]);

%figure(4);
%plot(fq, angle(S_nu(:,1,2)), fq, angle(S_x2(:,1,2)), ufq, angle(S_us(:,1,2)), '-o');
%legend('nufft-bia','nufft-x2','unif');  xlim([0,fq_max]);

%figure(5);
%plot(fq, abs(S_x2_bk(:,1,2)), '-+', fq, abs(S_x2(:,1,2)), '-o');

%figure(6);
%plot(fq, abs(S_x2_bk(:,1,1))-abs(S_x2(:,1,1)), '-o');


nGrangerT(uX,30)
tic();
getGCSapp(S_us)
tocs('GCapp us');
tic();
getGCSapp(S_nu)
tocs('GCapp nu');
tic();
getGCSapp(S_x2)
tocs('GCapp x2');

%S_mt = mX2S_mt(mX);
%SGrangerS(S_us)
%SGrangerS(S_mt)

%SGrangerS(S1)

%S2 = S1;
%S2(:,1,2) = S_us(:,1,2);
%S2(:,2,1) = S_us(:,2,1);
%getGCSapp(S2)
%SGrangerS(S2)

%getGCSapp(S1)



%df = @(S) S(2:end/2,1,1) - flipud(S(end/2+2:end,1,1));
%df21 = @(S) S(2:end/2,2,1) - flipud(S(end/2+2:end,1,2));
%df22 = @(S) S(2:end/2,2,1) - conj(S(2:end/2,1,2));

%%plot(S_us(2:end/2,1,1));
%%hold on
%%plot(flipud(S_us(end/2+2:end,1,1)), 'r');
%%hold off

%d1 = df(S_x2);
%plot(real(d1))

%plot(abs(df21(S1)))
%hold on
%plot(abs(df22(S1)), 'r')
%hold off

%plot(S1(:,1,1));


%fq_cut = 0.5;
%[S3, fqs3] = FreqCut(S1,fqs,fq_cut);
%[S_us3, fqs3] = FreqCut(S_us,fqs,fq_cut);
%[S_nu3, fqs3] = FreqCut(S_nu,fqs,fq_cut);

%figure(7);
%plot(fftshift(fqs), fftshift(S1(:,1,1)));

%figure(8);
%plot(fftshift(fqs3), fftshift(S3(:,1,1)));
%xlim([-fq_cut,fq_cut]);

%load('Ss.mat');

