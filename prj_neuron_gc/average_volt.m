% average of each group

%example:
% output_name = 'data/volt_net_3g_06_sc=0.01_pr=1_ps=0.007_t=5.00e+07_stv=0.5.dat';
% id = [1,1,1,2,2,2,3,3,3];  % the numbers are group id starting from 1, 0 means not to use.
% X = average_volt(output_name, id)
% save('-v7.3', 'X_3g_06.mat', 'X');  % use v7.3 so we can save data which is greater than to 2 GB on 64-bit systems in matlab.

function Y = average_volt(output_name, id)

% read from file or already the data itself
if ischar(output_name)
  p = length(id);
  fid = fopen(output_name, 'r');
  X = fread(fid, [p, Inf], 'double');
  fclose(fid);
else
  X = output_name;
end
len = size(X,2);

pp = max(id);
Y = zeros(pp, len);
for jj=1:pp
    Y(jj,:) = mean(X(id==jj,:),1);
end
clear('X');

end
