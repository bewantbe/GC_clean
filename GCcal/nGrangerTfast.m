% Multi-variate conditional Granger causality calculation in time domain.
% Much faster than pos_nGrangerT2() for large (>50) variable case.
% Almost as stable as pos_nGrangerT2(), and mathematically equivalent to it.
%
% Usage (see also nGrangerT()):
%   [GC, D, A2d] = nGrangerTfast(X, m, b_whiten_first)
%
% Time cost: O( len*m*p^2 ) + O( (p*m)^3 )
% RAM  cost: O( len*p ) + O( 3.5*(p*m)^2 ) * 8 Byte

function [GC, D, A2d] = nGrangerTfast(X, m, b_whiten_first)
  [p, len] = size(X);
  if p*p*m > p*(len-m)
    if len < p
      disp('Do you mis-transpose the data matrix?');
    end
    warning('Data length too short for GC.');
    GC = zeros(p, p);
    D = GC;
    A2d = [];
    return
  end
  if exist('b_whiten_first', 'var') && b_whiten_first ~= 0
    % whiten data can help reduce condition number
    switch b_whiten_first
    case 1
      % whiten time series data directly
      X = WhiteningFilter(X, m);  % whiten data, could change actual order needed.
      covz = getcovzpd(X, m);
    case 2
      % whiten the covariance in frequency domain
      fftlen = max(1024, 5*m);
      covz = getcovzpd(X, m);
      [A2d, D] = ARregressionpd(covz, p);
      S = A2S(A2d, D, fftlen);
      S = StdWhiteS(S);
      R = S2cov(S, m);
      covz = R2covz(R);
    otherwise
      error('no this mode');
    end
  else
    covz = getcovzpd(X, m);
  end

  [GC, D, A2d] = RGrangerTfast(covz, p);
end
