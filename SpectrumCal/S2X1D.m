% Spectrum Factorization of One Dimensional Case

function X = S2X1D(S)
  if (min(S)<=0) || iscomplex(S)
    warning('Error: can not factorize such spectrum!');
    X = NaN(size(S));
    return;
  end
  X  = exp(hilbert(0.5*log(S)));
end