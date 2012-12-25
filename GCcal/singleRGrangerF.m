% Calculate single directional group Granger Causality in frequency domain (by LSM)
% Possible to set different order for auto- and joint-regression
% very very slow version
%{
% Example:
    X = gdata(1e5, 3, 5);                  % get data: len=1e5, 3-od, 5-var
    R = getcovpd(X, 3);                    % get covariances up to order 3
    wGc = singleRGrangerF(R, [1 3], [2]);  % Causality from the [1 3]-th to the 2-th variable
    sprintf('GC=%19.16f', mean(wGc))       % check the result
%}
% 2012-01-21 xyy
% 2012-03-15 add od stuff; move all the actual calculation to singleGrangerFA().

function wGc = singleRGrangerF(R, id_y, id_x, od, fftlen)
if (exist('fftlen','var')==0)
  fftlen = 1024;
end

p=size(R,1);
lx = length(id_x);  ly = length(id_y);  lz = p-lx-ly;

if ~exist('od','var') || isempty(od)
  od_joint = round(size(R,2)/p)-1;
  od_auto  = od_joint;
elseif length(od)==1
  od_joint = od;
  od_auto  = od_joint;
elseif length(od)==2
  od_joint = od(1);
  od_auto  = od(2);
else
  error('unknown parameter type: od');
end

idx0 = true(1,p);  idx0(id_y)=false;
idx = true(1,(1+od_auto)*p);
idx(idx0(rem((1:((1+od_auto)*p))-1,p)+1)==false) = false;   % exclude y
[D, De0_star] = ARregression(R(idx0,idx));             % autoregression(without y)
[B, De0] = ARregression(R(:,1:p+p*od_joint));          % jonit-regression

wGc = singleGrangerFA(D, De0_star, B, De0, id_y, id_x, fftlen);

end

