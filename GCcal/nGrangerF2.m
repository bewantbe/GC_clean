% Granger Causality in frequency domain. only a prototype
% Example:
%{
    X = gdata(1e5, 2, 2);           % get data: len=1e5, 2-od, 2-var
    od = 6;
    [gcy2x, gcx2y, fx2y, fy2x, Sxx, Syy] = nGrangerF2(X, od);
    GC  = nGrangerT(X, od)          % check the results
    mean([fy2x fx2y])               % they must exactly the same(mathematically)
    plot(fy2x);
%}

function [gcy2x, gcx2y, fx2y, fy2x, Sxx, Syy] = nGrangerF2(X, od, fftlen)
if (exist('fftlen', 'var')==0)
	fftlen = 1024;
end
if (exist('od', 'var')==0)
	od = 64;
end

[gc,de,A] = pos_nGrangerT2(X, od);
%[A, noisecov] = ARregression(getcovpd(X, od));

[gcy2x, gcx2y, fx2y, fy2x, Sxx, Syy] = nGrangerFA2(A, de, fftlen);

end
