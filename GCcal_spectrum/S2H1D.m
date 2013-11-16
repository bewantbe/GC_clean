% Spectrum Factorization of One Dimensional Case

function [H, de] = S2H1D(S)
  [X, de] = S2X1D(S);
  H  = X/sqrt(de);
end
