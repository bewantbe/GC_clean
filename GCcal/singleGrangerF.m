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

%% solve in time domain
idx0 = true(1,p);  idx0(id_y) = false;
[~, De0_star, D] = pos_nGrangerT2(X(idx0,:), od);    % autoregression(without y)
[~, De0,      B] = pos_nGrangerT2(X, od);            % jonit-regression

wGc = singleGrangerFA(D, De0_star, B, De0, id_y, id_x, fftlen);

end

