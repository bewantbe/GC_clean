%
if strcmp(OCTAVE_VERSION,'3.4.3')==1
  addpath('/home/xyy/matlabcode/NUFFT_code/NUFFT_code/old_correct/NUFFT_code');
else
  addpath('/home/xyy/matlabcode/NUFFT_code/NUFFT_code/NUFFT_code');
end
tocs = @(st) fprintf('%s: t = %6.3fs\n', st, toc());

len = 1e5;
X = gdata(len,2,2);

mlen = 512;                                % length of each trial
n_trials = floor(size(X,2)/mlen);          % number of trials
mX = reshape(X(:,1:mlen*n_trials), size(X,1), mlen, []);
mT = fftshift(-mlen/2:mlen/2-1)'*ones(1,n_trials);

mX2 = mX;
mT2 = mT;

ufq = (0:mlen-1)/mlen;

tic; S_ft = mX2S_ft(mX);   tocs('ft');
tic; S_mt = mX2S_mt(mX);   tocs('mt');
tic; S_nuft = mX2S_nuft(mX, mT);  tocs('nuft');
tic; S_nuft_x2 = mX2S_nuft_x2(mX, mX2, mT, mT2);  tocs('nuft_x2');

mad = @(xx) mean(abs(xx(:)));

mad(S_ft-S_mt)
mad(S_ft-S_nuft)
mad(S_ft-S_nuft_x2)

figure(1);
plot(ufq, S_ft(:,1,1));
figure(2);
plot(ufq, S_mt(:,1,1));

