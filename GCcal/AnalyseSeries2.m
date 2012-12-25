% Analyse the data in detail
%e.g. (see analyse_GC_simple.m for more detailed example)
%{
    addpath ..
    len = 1e5;
    X = gdata(len, 2, 3);
    s_od = 1:99;
    [oGC, oDe, R] = AnalyseSeries(X, s_od);
    [aic_od, bic_od, zero_GC, oAIC, oBIC] = AnalyseSeries2(s_od, oGC, oDe, len);
    disp('GC(od=20):');
    disp(oGC(:,:,20));
%}

function [aic_od, bic_od, zero_GC, oAIC, oBIC] = AnalyseSeries2(s_od, oGC, oDe, len)

p = size(oGC, 1);
oAIC = zeros(1, length(s_od));
oBIC = zeros(1, length(s_od));
for k=1:length(s_od)
    npm = p^2 * s_od(k);          % number of free parameters
    oAIC(k) = log(det(oDe(:,:,k))) + 2 * npm / len;
    oBIC(k) = log(det(oDe(:,:,k))) + npm * log(len) / len;
%??    oAICc(k)= log(det(oDe(:,:,k))) + 2 * npm / len + 2*npm*(npm+1)/(len-npm-1)/len;
end

[~, od_id]  = min(oBIC);
bic_od = s_od(od_id);
[~, od_id]  = min(oAIC);
aic_od = s_od(od_id);

[~, od_id] = min(abs(s_od-min(aic_od, max(s_od)-10)));    % at least 10 average
% linear interpretation to find out GC without order bias
mid_GC = mean(oGC(:,:,od_id:end), 3);   % middle value
mid_od = mean(s_od(od_id:end));
%zero_GC = abs(mid_GC - mid_od/len);
zero_GC = mid_GC - mid_od/len;
zero_GC(1:p+1:p*p) = 0;

end
