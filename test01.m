%

%X = gdata(1e6,2,2);

%n_trials = 100;
%mlen = floor(size(X,2)/n_trials);         % length of each trials at most
mlen = 512;
n_trials = floor(size(X,2)/mlen);          % number of trials at most
tic;
mX = reshape(X(:,1:mlen*n_trials), size(X,1), mlen, []);
toc

%tic;
%S = mX2S(mX);
%toc

%tic;
%S = mX2S_ft(mX);
%toc

%plot(S(:,1,1));
