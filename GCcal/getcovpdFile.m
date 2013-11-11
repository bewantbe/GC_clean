% Calculate covariance matrices
% Prucde the same result as getcovpd()
% But cost less memory

% X is a multivariate time series p*len dimension
% p is the number of variates, len is the length of time
% m is the maximun offset of covariance to calculate
% result R is a p * (p*(m+1)) dimension matrix

% you may specify skip how many data points(each 8*p byte) by x_xkip

function R = getcovpdFile(XF, p, m, x_skip, x_len)

if (~exist('XF','var') || ~exist('p','var')|| ~exist('m','var') || isempty(XF))
  disp('usage: R = getcovpdFile(XFileName, p, m)');
  error('XFileName must be file name to the X matrix');
end
if ~ischar(XF)
  error('XFileName should be a string (file name)');
end
if ~exist(XF, 'file')
  error('File not found.');
end
d = stat(XF);
f_len = round(d.size/8/p);
if abs(f_len*p*8-d.size)>0.5
  error('Data file size does not compatible with len and p');
end
if ~exist('x_skip', 'var')
  x_skip = 0;              % read from the beginning
end
if ~exist('x_len', 'var')
  len = f_len - x_skip;    % read to the end
end
if len<0 || x_skip<0
  error('Length or skip should be positive number');
end
if x_skip+len>f_len
  error('Request data too many');
end
if len<=m || m<0
  error('No enough data to calculate such order');
end

fid = fopen(XF, 'r');
if x_skip>0
  fseek(fid, x_skip*p*8);
end

mem_lim = 2^29;            % memory limit in byte
len_seg = round(mem_lim/8/p);

% we need to make sure that X_head and X_tail are non-overlapped
% hence if the length is short enough, directly calculate it by getcovpd
if ceil(len/len_seg)==1
  X = fread(fid, [p, len_seg], 'double');
  fclose(fid);
  R = getcovpd(X, m);
  return
end
X_head = fread(fid, [p, m], 'double');  % in ordre match original result we need this part
X      = X_head;
X_glue = zeros(p, m);
R      = zeros(p, p*(m+1));
aveX   = sum(X,2);

for k=1:ceil(len/len_seg)
  X_glue = X(:, end-m+1:end);
  clear('X');                   % hope to help save memory...
  X = fread(fid, [p, len_seg], 'double');
  aveX = aveX + sum(X,2);
  X = [X_glue, X];
  for k = 0 : m
    R(:, 1+p*k : p+p*k) = R(:, 1+p*k : p+p*k) + (X(:, m+1:end) * X(:, m-k+1:end-k)');
  end
end
fclose(fid);

R = R / (len-m);
aveX = aveX / len;
for k = 0 : m
  R(:, 1+p*k : p+p*k) = R(:, 1+p*k : p+p*k) - aveX*aveX';
end

% make the result match getcovpd exactly (mathematically)
X_tail = X(:, end-m+1:end);
R_comp = zeros(p, p*(m+1));
for k = 0 : m
  R_comp(:, 1+p*k : p+p*k) = (sum(X_head,2) - m*aveX)*aveX'...
    + aveX*(sum([X_head(:,1:m-k), X_tail(:,end-k+1:end)],2) - m*aveX)';
end
R = R + R_comp/(len-m);
