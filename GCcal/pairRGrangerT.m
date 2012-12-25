% Compute pair-wise Granger Causality of p-variable input

% time cost is about: O(p^2 *8*m^2)

function GC = pairRGrangerT(R)
[p,l] = size(R);

GC = zeros(p,p);
for ii=1:p
    for jj=ii+1:p
	idx = [ii:p:l; jj:p:l];
        gc = RGrangerT([R(ii,idx); R(jj,idx)]);
        GC(jj,ii) = gc(2,1);
        GC(ii,jj) = gc(1,2);
    end
end

end

