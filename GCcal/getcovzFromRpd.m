%
% Mean value of X must be subtracted.

function covz = getcovzFromRpd(X, R, m)
[p, len] = size(X);
if p~=size(R,1)
  error('wrong paire of X R');
end
m_R = round(size(R,2)/p)-1;
if m_R < m
  error('Order of R is too low');
end
% Here assume R is obtained by
% R = getcovpd(X, m_R);

t_bg = m+1;     % time range in common
t_ed = len-m;
% middle part
R  = (len-m_R) * R;
RR = zeros(p, (2*m+1)*p);
for k = 0 : m
  %R_k_ans = X(:, t_bg:t_ed) * X(:, t_bg-k:t_ed-k)';
    R_hd = X(:, t_bg:m_R) * X(:, t_bg-k:m_R-k)';
    R_ed = X(:, t_ed+1:len) * X(:, t_ed+1-k:len-k)';
    R_k = R(:,k*p+1:(k+1)*p) + (R_hd - R_ed);
  %max(abs((R_k_ans-R_k))(:))
    RR(:, (m+k)*p+1 : (m+1+k)*p) = R_k;
    RR(:, (m-k)*p+1 : (m+1-k)*p) = R_k';
end
% head and tail parts
covz = getcovzpdhded(X,t_bg,t_ed,m);
% combine them
for k=0:m
    covz(k*p+1:(k+1)*p,:) = covz(k*p+1:(k+1)*p,:) + RR(:,(m-k)*p+1:(2*m+1-k)*p);
end
covz = covz/(len-m);
