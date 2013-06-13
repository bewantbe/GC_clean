%

%X = gdata(1e6,2,2);

X = gendata_neu('net_2_2', 0.02, 1, 0.012, 1e6, 0.5);
X = bsxfun(@minus, X, mean(X,2));

mlen = 512;                                % length of each trials at most
n_trials = floor(size(X,2)/mlen);          % number of trials at most
mX = reshape(X(:,1:mlen*n_trials), size(X,1), mlen, []);

S_a = mX2S_ft(mX);
S_mt = mX2S_mt(mX);

S = permute(S_a,[3,1,2]);
figure(1);
semilogy(S(:,1,1));
figure(2);
semilogy(abs(S(:,1,2)));

tic
gc = nGrangerT(X,20)
toc

tic
gc_np_a = SGrangerS(S_a)
toc

tic
gc_np_mt = SGrangerS(S_mt)
toc

