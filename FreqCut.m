% lowpass cut

function [S3 fqs3] = FreqCut(S1, fqs, fq_cut)

  %ids = fq_cut < fq & fq < (fq_max-fq_cut);
  slen = size(S1,1);
  fq_max = abs(fqs(end/2));
  l1 = floor(slen/2*fq_cut/fq_max);
  ids = l1+1 : slen-l1;
  S3 = ipermute(S1, [3,1,2]);   % convert to p*p*fftlen
  S3(:,:,ids) = [];
  S3 = permute(S3, [3,1,2]);   % convert to fftlen*p*p
  fqs3 = fqs;
  fqs3(ids) = [];

end
