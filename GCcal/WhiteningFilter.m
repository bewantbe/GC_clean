%% Get residuals and coefficients in auto-regression of each variables.

function [srd, sas] = WhiteningFilter(X, use_od, coeff_mode)
[p, len] = size(X);
if exist('coeff_mode','var') && strcmpi(coeff_mode,'qr')
    FitCoef = @(x) pos_nGrangerT_qr(x, use_od);
else
    FitCoef = @(x) pos_nGrangerT2(x, use_od);
end

srd = zeros(p, len);
sas = zeros(p, use_od);
for k=1:p
    [~, ~, sa] = FitCoef(X(k,:));
    srd(k,:) = filter([1, sa], [1], X(k,:)-mean(X(k,:),2));
    sas(k,:) = sa;
end
%srd(:,1:use_od) = 0;     % there is no data for this part

end
