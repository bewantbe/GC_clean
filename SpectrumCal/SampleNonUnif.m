% Obtain non-uniform samples from uniform samples

%Input
% uX : p * len
% mlen: (original) sample length of each trial (sample window size)
% slen: length of resample data

%Output
% mX1: p * slen * n_trials
% mT1: slen * n_trials

function [mX1, mX2, mT1, mT2] = SampleNonUnif(uX, mlen, slen, smode, more_trials)
  n_trials = floor(size(uX,2)/mlen);   % number of nonoverlap trials at most
  uX = bsxfun(@minus, uX, mean(uX,2));

  p = size(uX,1);
  if ~exist('smode','var')
    smode='1';
  end
  if isnumeric(smode)
    more_trials = smode;
    smode='1';
  end
  l_inc = mlen;    % "sample" window increment
  if exist('more_trials','var')
    n_trials = more_trials;
    if (n_trials>1)
      l_inc = (size(uX,2)-mlen)/(n_trials-1);  % may not an integer
    end
  end
  switch smode
    case '0'  % (approximate) equal space
      mT1 = repmat(round(1 : mlen/slen : mlen)', [1, n_trials]);
      mT2 = mT1;
    case '1'  % uniformly random
      % sort is not necessary but safe to add it
      mT1 = sort(ceil(mlen*rand(slen,n_trials)),1);  
      mT2 = sort(ceil(mlen*rand(slen,n_trials)),1);
    case '2'  % equal space with random perturbation
      dt  = mlen/slen;
      T   = repmat((1:dt:mlen)', [1, n_trials]);
      mT1 = floor(T + dt*rand(slen, n_trials));
      mT2 = floor(T + dt*rand(slen, n_trials));
  end
  mX1 = zeros(p,slen,n_trials);
  mX2 = zeros(p,slen,n_trials);
  for tr=1:n_trials
    mX1(:,:,tr) = uX(:, mT1(:,tr)+round(l_inc*(tr-1)));
    mX2(:,:,tr) = uX(:, mT2(:,tr)+round(l_inc*(tr-1)));
  end
end
