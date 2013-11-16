% Calculate Granger Causality through spectrum

function [GC, de, H] = nGrangerS(X, mlen)
if ~exist('mlen','var')
    mlen = 512;                            % length of each trials at most
end

X = bsxfun(@minus, X, mean(X,2));

n_trials = floor(size(X,2)/mlen);          % number of trials at most
mX = reshape(X(:,1:mlen*n_trials), size(X,1), mlen, []);

S_mt = mX2S_mt(mX);

[GC, de, H] = SGrangerS(S_mt);

end
