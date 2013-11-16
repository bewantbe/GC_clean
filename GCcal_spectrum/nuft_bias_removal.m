% Remove the bias caused by non-uniform Fourier Transform in Spectrum

function S2 = nuft_bias_removal(S, De, stv, slen)
  p = size(S,2);
  for j=1:p
    for k=1:p
      S2(:,j,k) = (S(:,j,k) - De(j,k)) * stv * slen/(slen-1);
    end
  end
  %S1(:,1,1) = S1(:,1,1) - De(1,1)*stv;
  %S1(:,2,2) = S1(:,2,2) - De(2,2)*stv;
  %S1(:,1,2) = S1(:,1,2) - De(1,2)*stv;
  %S1(:,2,1) = S1(:,2,1) - De(2,1)*stv;
end
