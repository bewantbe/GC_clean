% Compute pair-wise Granger Causality of p-variable input

% time cost is about: O(len * m * p^2 + p^2 *8*m^2)

function GC = pairGrangerT(X, m, bad_mode)
if (nargin ~= 2 && nargin ~= 3)
    error('Usage: pairGrangerT(X, m, bad_mode), m is the order of AR');
end
if ~exist('bad_mode','var')
  bad_mode = 0;
end

switch bad_mode
case 0
  GC = pairRGrangerT(getcovpd(X, m));
case 1
  p = size(X, 1);
  GC = zeros(p,p);
  for ii=1:p
      for jj=ii+1:p
          gc = pos_nGrangerT_qrm(X([ii,jj],:), m);
          GC(jj,ii) = gc(2,1);
          GC(ii,jj) = gc(1,2);
      end
  end
end

