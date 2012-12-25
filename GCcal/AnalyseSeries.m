% analysis a given time series in different order
% input:
%   X: p*len matrix
%   s_od: orders to calculate, e.g. s_od = 1:99
% For very bad condition problem you may try this:
%   [oGC, oDe, R] = AnalyseSeries(X, s_od, 1)

function [oGC, oDe, R] = AnalyseSeries(X, s_od, bad_mode)
if ~exist('bad_mode','var')
  bad_mode = 0;
end

p = size(X, 1);

% calculate GC in different order
R = getcovpd(X, max(s_od));
oGC = zeros(p, p, length(s_od));
oDe = zeros(p, p, length(s_od));
for k=1:length(s_od)
  if ~bad_mode
    [cGC Deps] = RGrangerT(R(:, 1:p+p*s_od(k)));
  else
    [cGC Deps] = pos_nGrangerT2(X, s_od(k));   % for bad condition problem
  end
  oGC(:, :, k) = cGC;
  oDe(:, :, k) = Deps;
end

end
