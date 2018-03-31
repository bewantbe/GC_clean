% solve partial GC problem
%  use y only to help predict x
%  return log( var(x) / var(x-Axdy*y) )

% e.g.
%  e1=randn(1,1e5);
%  e2=randn(1,1e5);
%  partial_GC(e1(2:end)+e2(1:end-1), e2(2:end), 2)
% should give result close to log(2)

function [pGC, jrdx] = partial_GC(srdx, srdy, od)
len = length(srdx);
if (exist('od', 'var') == 0)
  od = 20;
end

covxy = zeros(1,od);
for ii=1:od
    covxy(ii) = (srdx(od+1:end) * srdy(od-ii+1:end-ii)') / (len-od);
end

%{
Z = zeros(od,len-od);
for l = 1 : od
    Z(l,:) = srdy(od-l+1:len-l);
end
covyy = Z*Z'/(len-od);
%}
covy = getcovpd(srdy, od-1);
long_covy = [fliplr(covy(2:end)),covy];
covyy = zeros(od, od);
for ii=1:od
    covyy(ii,:) = long_covy(od-ii+1:end-ii+1);
end
Axdy = (covyy \ (-covxy'))';

jrdx = srdx + filter([0, Axdy], [1], srdy);

pGC = log( var(srdx) / var(jrdx) );

end
