% Granger Causality in frequency domain (by LSM).
% a very general but slow procedure.
% returen value wGC is a p*p*fftlen matrix
%{
% Example:
    X = gdata(1e5, 3, 5);           % get data: len=1e5, 3-od, 5-var
    [od_joint,od_auto]=chooseOrderFull(X);
    od = max([od_joint, od_auto])   % choose a high enough order
    wGC = nGrangerF(X, od, 1024);   % frequency domain GC, use fftlen=1204
    GC  = nGrangerT(X, od)          % check the results
    mean(wGC,3)                     % they must exactly the same(mathematically)
    plot(squeeze(wGC(2,1,:)));
%}
% 2012-01-21 xyy

function wGC = nGrangerF(X, od, fftlen, bad_mode)
if (exist('fftlen', 'var')==0)
	fftlen = 1024;
end
if (exist('od', 'var')==0)
	od = 64;
end
if ~exist('bad_mode','var')
	bad_mode = 0;
end
if length(od)>1
  od_joint = od(1);  % order for joint AR fitting
  od_R = max(od);      % order for cov series
else
  od_joint = od;
  od_R = od;
end

if bad_mode > 0
  [gc, De, As2] = pos_nGrangerT_qrm(X, od_joint);
  S = A2S(As2(:,1:end-1), De, fftlen);
  R = S2cov(S, od_R);
else
  R = getcovpd(X, od_R);
end
p = size(R,1);
wGC = zeros(p,p,fftlen);
for ki=1:p
  for kj=1:p
    if ki==kj
      continue;
    end
    wgc = singleRGrangerF(R, kj, ki, od, fftlen);
    wGC(ki,kj,:) = wgc;
  end
end

end
