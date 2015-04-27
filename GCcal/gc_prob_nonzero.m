% Calculate the non-zero probability of each GC. Also used in conditional GC.
% only suitable for nGrangerT series functions(use same fitting order
%  in auto and joint regression).
%
% For case of small sample, the following give more accurate result
% nzp = gc_prob_nonzero(sv, od, len-(p+1)*od);

function p = gc_prob_nonzero(gc, od, len)
if nargin<3
  disp(' Calculate the non-zero probability of each GC. Can also used in conditional GC.');
  usage('p = gc_nonzero_prob(gc, od, len)');
end
p = chi2cdf(gc*len, od);
end
