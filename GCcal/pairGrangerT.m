% Pairwise Granger causality calculation in time domain.
%
% Time cost is about: O(len * m * p^2) + O(p^2 *16*m^3)

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
  return
case 1
  f_gc = @(e_X, e_m) pos_nGrangerT2(e_X, e_m);
case 2
  f_gc = @(e_X, e_m) pos_nGrangerT_qrm(e_X, e_m);
end

p = size(X, 1);
GC = zeros(p,p);
for ii=1:p
    for jj=ii+1:p
        gc = f_gc(X([ii,jj],:), m);
        GC(jj,ii) = gc(2,1);
        GC(ii,jj) = gc(1,2);
    end
end

