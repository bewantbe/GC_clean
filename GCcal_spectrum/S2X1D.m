% Spectrum Factorization of One Dimensional Case

function [X, de] = S2X1D(S)
  if ndims(S)>2
    error('S2X1D: Only accept 1-Dimensional data');
  end
  if (min(S)<=0) % || iscomplex(S), no iscomplex in matlab
    warning('Error: can not factorize such spectrum! (shold have S>0)');
    X = NaN(size(S));
    de= NaN;
    return;
  end
  X  = exp(conj(hilbert(0.5*log(S))));
  de = abs(mean(X))^2;
end
