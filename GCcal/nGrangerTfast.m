%

function GC = nGrangerTfast(X, m, b_whiten_first)
  p = size(X, 1);
  if exist('b_whiten_first', 'var') && b_whiten_first ~= 0
    % whiten data can help reduce condition number
    switch b_whiten_first
    case 1
      % whiten time series data directly
      X = WhiteningFilter(X, m);  % whiten data
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

  GC = RGrangerTfast(covz, p);
end
