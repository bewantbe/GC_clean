% Pairwise Granger causality calculation in time domain.
%
% Time cost is about: O(p^2 *16*m^3)

function [pairGC, de, pairA] = pairRGrangerT(R)
[p,l] = size(R);
m = round(l/p)-1;

pairGC = zeros(p,p);
pairA = zeros(p, p*m);
for ii=1:p
    for jj=ii+1:p
	idx = [ii:p:l; jj:p:l];
        [gc, ~, a] = RGrangerT([R(ii,idx); R(jj,idx)]);
        pairGC(ii,jj) = gc(1,2);
        pairGC(jj,ii) = gc(2,1);
        pairA(ii,jj:p:end) = a(1,2:2:end);
        pairA(jj,ii:p:end) = a(2,1:2:end);
    end
end

de = nan;

end

