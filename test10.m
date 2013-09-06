%

len = 1e5;
len_ext = 1000;
noisecov = [1.0, 0.4; 0.4, 0.7];
A = [-0.9 ,  0.0, 0.5, 0.0;
     -0.16, -0.8, 0.2, 0.5];
noisecov = diag([0.3 1.0 0.2]);
A = [-0.8  0.0 -0.4  0.5 -0.2  0.0;
      0.0 -0.9  0.0  0.0  0.8  0.0;
      0.0 -0.5 -0.5  0.0  0.0  0.2];

tic; X=gendata_linear(A,noisecov,1e6); toc; figure(2); hist(X(1,:),-9:0.1:9)
tic; X=gendata_linear(A,noisecov,1e7); toc; tic; [gc,ede,eA]=nGrangerT(X,2); toc; eA

%oX = gendata_linear(A,noisecov,len+len_ext);
%oX = oX(:,len_ext+1:end);

oX = gdata(len,2,2);

stv0 = 1;
slen = 128;                           % points in one non-unif sample

T_segment = stv0*slen;

stv = stv0;
uX = oX(:,1:(stv/stv0):end);
mlen = T_segment/stv;                 % length of each trials
n_trials = floor(size(uX,2)/mlen);    % number of trials at most
mX = reshape(uX(:,1:mlen*n_trials), size(uX,1), mlen, []);
ufq = (0:mlen-1)/mlen/stv;
fq_max = slen/T_segment;

S_us = mX2S_ft(mX);


figure(2);
plot(ufq, real(S_us(:,1,1)), '-o');
xlim([0,fq_max]);

