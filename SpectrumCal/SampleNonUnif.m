% Obtain non-uniform samples from uniformly sampled data

%Input
% uX : p * len
% mlen: (original) sample length of each trial (sample window size)
% slen: length of resample data

%Output
% mX: p * slen * n_trials
% mT: slen * n_trials

function [mX, mT] = SampleNonUnif(uX, mlen, slen, smode, more_trials)
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
    case {'0' 'u'}  % (approximate) equal space
      mT = repmat(round(1 : mlen/slen : mlen)', [1, n_trials]);
    case {'1' 'r'}  % uniformly random
      % sort is not necessary but safe to add it
      mT = sort(ceil(mlen*rand(slen,n_trials)), 1);  
    case '2'  % equal space with random perturbation
      dt  = mlen/slen;
      T   = repmat((1:dt:mlen)', [1, n_trials]);
      mT = floor(T + dt*rand(slen, n_trials));
    otherwise
      error("Unsupported resample mode");
  end
  mX = zeros(p,slen,n_trials);
  for tr=1:n_trials
    mX(:,:,tr) = uX(:, mT(:,tr)+round(l_inc*(tr-1)));
  end
end
