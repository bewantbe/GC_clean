%% get residual and coefficients of self-regression

function [srd, sas] = WhiteningFilter(X, use_od)
[p, len] = size(X);

srd = zeros(p, len);
sas = zeros(p, use_od);
for k=1:p
    [~, ~, sa] = pos_nGrangerT2(X(k,:), use_od);
    srd(k,:) = filter([1, sa], [1], X(k,:)-mean(X(k,:),2));
    sas(k,:) = sa;
end
srd(:,1:use_od) = 0;     % there is no data for this part

end
