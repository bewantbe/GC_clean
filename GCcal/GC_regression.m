% General GC regression analysis on a single realization (data)

% [srd, rd, aveX] = GC_regression(X[, od_max[, od]])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input data
[p,len] = size(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis
if ~exist('od_max', 'var')
  od_max = 99;
end
s_od = 1:od_max;
%[oGC, oDe, R] = AnalyseSeries(X, s_od);
[oGC, oDe, R] = AnalyseSeriesFast(X, s_od);
[aic_od, bic_od, zero_GC, oAIC, oBIC] = AnalyseSeries2(s_od, oGC, oDe, len);
%save('-v6','tmp_gc.mat','p','len','s_od','oGC','oDe','R');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get residual (p-var)
aveX = mean(X, 2);
X = bsxfun(@minus, X, aveX);

od = aic_od;                               % order used to get residuals
%  [oa, os] = chooseOrderFull(X, 'AIC', 99);
%  od = max([oa, os]);
srd = zeros(p,len);
for k=1:p
  sa = ARregression(R(k,k:p:p+p*od));
  srd(k,:) = filter([1, sa], [1], X(k,:));
end

A  = ARregression(R(:,1:p+p*od));
rd = zeros(p,len);
if p==1
  for k=1:p
    rd(k,:) = filter([1, A(k,k:p:end)], [1], X(k,:));
    for j=1:p
      if (k==j)
        continue;
      end
      rd(k,:) = filter([0, A(k,j:p:end)], [1], X(j,:)) + rd(k,:);
    end
  end
else
  rd = MAfilter_v5(A, X);
end

X = bsxfun(@plus, X, aveX);

