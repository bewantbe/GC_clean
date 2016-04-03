% Calculate single directional group Granger Causality in frequency domain (by LSM)
% very very slow version
%{
% Example:
    X = gdata(1e5, 3, 5);                  % get data: len=1e5, 3-od, 5-var
    wGc = singleGrangerF(X, [1 3], [2]);  % Causality from the [1 3]-th to the 2-th variable
    sprintf('GC=%19.16f', mean(wGc))       % check the result
%}
% 2012-01-21 xyy
% 2012-03-15 xyy: move all the actual calculation to singleGrangerFA().

function wGc = singleGrangerF(X, id_y, id_x, od, fftlen)
if ~exist('fftlen','var')
  fftlen = 1024;
end
if ~exist('od','var')
  od = 64;
end

[p, len]=size(X);

tmp_v = zeros(1,p);  tmp_v(id_y) = 1;                % can not overlap between groups
if sum(tmp_v(id_x))>0 || length(id_x)+length(id_y)>=p
  disp('id_x:');  disp(id_x);
  disp('id_y:');  disp(id_y);
  error('repeated index!!');
end

% Solve in time domain
covz = getcovzpd(X, od);
id_no_y = true(p, 1);  id_no_y(id_y) = false;
id_no_y = repmat(id_no_y, 1, od+1);
[D, De0_star] = ARregressionpd(covz(id_no_y, id_no_y), p-length(id_y));
[B, De0     ] = ARregressionpd(covz, p);
wGc = singleGrangerFA(D, De0_star, B, De0, id_y, id_x, fftlen);

end

