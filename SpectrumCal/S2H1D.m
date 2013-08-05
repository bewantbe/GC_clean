% Spectrum Factorization of One Dimensional Case

function [H, de] = S2H1D(S)
  if (min(S)<=0) || iscomplex(S)
    warning('Error: can not factorize such spectrum!');
    H = NaN(size(S));
    de= NaN;
    return;
  end
  X  = exp(hilbert(0.5*log(S)));
  de = abs(mean(X));  % or use real()?
  H  = X/de;
  de = de*de;
end
