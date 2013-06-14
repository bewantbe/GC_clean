%

%X = gdata(1e6,2,2);
X = gdata(1e6,2,3);

%X = gendata_neu('net_2_2', 0.02, 1, 0.012, 1e6, 0.5);
X = bsxfun(@minus, X, mean(X,2));

mlen = 512;                                % length of each trials at most
n_trials = floor(size(X,2)/mlen);          % number of trials at most
mX = reshape(X(:,1:mlen*n_trials), size(X,1), mlen, []);

S_a = mX2S_ft(mX);
S_mt = mX2S_mt(mX);

tic
gc = nGrangerT(X,20)
toc

tic
gc_np_a = SGrangerS(S_a)
toc

tic
gc_np_mt = SGrangerS(S_mt)
toc

S_a = permute(S_a,[3,1,2]);
S_mt = permute(S_mt,[3,1,2]);

fftlen = size(S_a,1);
fq = (0:fftlen-1)/fftlen;
cfq = fq - 0.5;

figure(1);
semilogy(cfq, fftshift(S_a(:,1,1)), cfq, fftshift(S_mt(:,1,1)));
xlim([-0.5,0.5]);
figure(2);
semilogy(cfq, fftshift(abs(S_a(:,1,2))), cfq, fftshift(abs(S_mt(:,1,2))));
xlim([-0.5,0.5]);

