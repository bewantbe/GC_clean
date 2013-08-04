% obtain non-uniform samples from uniform samples
% uX : p * len
% mX1: p * slen * n_trials
% mT1: slen * n_trials

function [mX1, mX2, mT1, mT2]=SampleNonUnif(uX, mlen, slen, smode)
  %mlen = 512;                            % length of each trial
  n_trials = floor(size(uX,2)/mlen);      % number of trials
  uX = bsxfun(@minus, uX, mean(uX,2));
  mX = reshape(uX(:,1:mlen*n_trials), size(uX,1), mlen, []);

  p = size(uX,1);
  %slen = 128;
  if ~exist('smode','var')
    smode='1';
  end
  switch smode
  case '1'
    % sort is not necessary but safe to add it
    mT1 = sort(round(mlen*rand(slen,n_trials)+0.5),1);  
    mT2 = sort(round(mlen*rand(slen,n_trials)+0.5),1);
  case '0'
    mT1 = repmat(round(1 : mlen/slen : mlen)', [1, n_trials]);
    mT2 = mT1;
  case '2'
    dt  = mlen/slen;
    T   = repmat((1:dt:mlen)', [1, n_trials]);
    mT1 = floor(T + dt*rand(slen, n_trials));
    mT2 = floor(T + dt*rand(slen, n_trials));
  end
  mX1 = zeros(p,slen,n_trials);
  mX2 = zeros(p,slen,n_trials);
  for tr=1:n_trials
    mX1(:,:,tr) = mX(:,mT1(:,tr),tr);
    mX2(:,:,tr) = mX(:,mT2(:,tr),tr);
    %mX1(:,:,tr) = uX(:,mT1(:,tr)+mlen*(tr-1));
    %mX2(:,:,tr) = uX(:,mT2(:,tr)+mlen*(tr-1));
  end
end
