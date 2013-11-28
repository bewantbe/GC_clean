% get positive defined COVZ

function covz = getcovzpd(X, m)

[p, len] = size(X);
% make sure that the average of X is zero
%X = X - mean(X,2)*ones(1,len);
X = bsxfun(@minus, X, mean(X,2));

t_bg = m+1;     % time range in common
t_ed = len-m;
% middle part
R = zeros(p, (2*m+1)*p);
for k = 0 : m
    R_k = X(:, t_bg:t_ed) * X(:, t_bg-k:t_ed-k)';
    R(:, (m+k)*p+1 : (m+1+k)*p) = R_k;
    R(:, (m-k)*p+1 : (m+1-k)*p) = R_k';
end
% head and tail parts
covz = zeros((m+1)*p,(m+1)*p);
for i1=0:m
    for i2=0:m
        k=i2-i1;
        if (k>=0)
            covz(i1*p+1:(i1+1)*p, i2*p+1:(i2+1)*p) =...
              X(:,m+1-i1:t_bg-1)*X(:,m+1-i1-k:t_bg-1-k)'...
            + X(:,t_ed+1:len-i1)*X(:,t_ed+1-k:len-i1-k)';
        else
            covz(i1*p+1:(i1+1)*p, i2*p+1:(i2+1)*p) =...
              X(:,m+1-i2+k:t_bg-1+k)*X(:,m+1-i2:t_bg-1)'...
            + X(:,t_ed+1+k:len-i2+k)*X(:,t_ed+1:len-i2)';
        end
    end
end
% combine them
for k=0:m
    covz(k*p+1:(k+1)*p,:) = covz(k*p+1:(k+1)*p,:) + R(:,(m-k)*p+1:(2*m+1-k)*p);
end
covz = covz/(len-m);
