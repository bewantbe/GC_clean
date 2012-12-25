% Choose Order
% ic_mode can be 'BIC', 'AIC' or 'AICc'

function [best_od, xic] = chooseOrder(X, ic_mode, od_max)
if ~exist('ic_mode','var')
    ic_mode = 'BIC';
end
if ~exist('od_max', 'var')
    od_max = 99;
end

R = getcovpd(X, od_max);
[best_od, xic] = chooseROrder(R, length(X), ic_mode);

end
