% Analyse a given time series in a series of orders
% Input & Output: See comments in AnalyseSeries.m
% Since this function relies on Levinson algorithm, it not suitable for
% bad condition problem (e.g. very short data, filtered data)
%
% Time cost (haven't verified):
%   O(p^4 * m^2) (from theory)

function [oGC, oDe, R] = AnalyseSeriesFast(X, s_od, b_input_cov)
  if exist('b_input_cov','var') && b_input_cov
    % so X is actually covariance
    R = X;
    p = size(R, 1);
    m = round(size(R, 2) / p) - 1;
    if (m < max(s_od))
      error('GC_clean: AnalyseSeriesLevinson: no enough data');
    end
  else
    p = size(X, 1);
    m = max(s_od);
    R = getcovpd(X, m);
  end
oGC = zeros(p, p, length(s_od));
oDe = zeros(p, p, length(s_od));

[~, s_De] = BlockLevinson(R);
for k=1:length(s_od)
  oDe(:,:,k) = s_De(:,:,s_od(k));
end

if p==1
  return;
end

dj = zeros(p,1);
for j=1:p   % omit j-th variable
  [~, s_Dej] = BlockLevinson(R((1:p)~=j, mod(0:p*(m+1)-1, p)+1~=j));
  % calculate GC from j-th variable in different order
  for k=1:length(s_od)
    dd = diag(oDe(:,:,k));
    dj((1:p)~=j) = diag(s_Dej(:,:,s_od(k)));
    dj(j) = dd(j);
    oGC(:, j, k) = log(dj ./ dd);
  end
end

end
