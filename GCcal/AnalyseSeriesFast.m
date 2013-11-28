% Analysis a given time series in different order
% Input:
%   X: p*len matrix
%   s_od: orders to calculate, e.g. s_od = 1:99
% Since this function relies on Levinson algorithm, it not suitable for
% vary bad condition problem

% Time cost (haven't verified):
%   O(m^1.5) (from numerical), or
%   O(p^4 * m^2) (from theory)

function [oGC, oDe, R] = AnalyseSeriesFast(X, s_od)
if ~exist('bad_mode','var')
  bad_mode = 0;
end
m = max(s_od);
R = getcovpd(X, m);
p = size(X, 1);
oGC = zeros(p, p, length(s_od));
oDe = zeros(p, p, length(s_od));

[~, s_De] = BlockLevinson(R);
for k=1:length(s_od)
  oDe(:,:,k) = s_De(:,:,s_od(k));
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
