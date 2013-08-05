% make mS Hermitian, and mS(-f) = mS(f).'

function mS = Makeup4SpectrumFact(S)
  S = ipermute(S, [3,1,2]);   % convert to p*p*fftlen
  len = size(S,3);
  p   = size(S,1);
  mS = zeros(p,p,len);
  mS(:,:,1) = 0.25*(S(:,:,1) + S(:,:,1).' + S(:,:,1)' + conj(S(:,:,1)));
  for k=2:(len/2+1)
    mS(:,:,k) = 0.25*( S(:,:,k        ) + S(:,:,k        )'
                    + (S(:,:,(len-k+2)) + S(:,:,(len-k+2))' ).' );
    mS(:,:,len-k+2) = mS(:,:,k).';
  end
  mS = permute(mS, [3,1,2]);
end
