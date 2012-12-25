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

function wGC = nGrangerF(X, od, fftlen)
if (exist('fftlen', 'var')==0)
	fftlen = 1024;
end
if (exist('od', 'var')==0)
	od = 64;
end

R = getcovpd(X, od);
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
