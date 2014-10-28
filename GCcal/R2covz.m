
function covz = R2covz(R)

[p, m] = size(R);
m = round(m/p);

% covariance series that has negative time
RR = zeros(p,p*(2*m-1));
RR(:, m*p-p+1:end) = R(:,1:m*p);
for k = 1 : m-1
    RR(:,(m-k)*p-p+1:(m-k)*p) = R(:,k*p+1:k*p+p)';
end
% construct the big covariance matrix
covz = zeros(p*m, p*m);
for k = 0 : m-1
    covz(k*p+1:k*p+p, :) = RR(:, (m-k)*p-p+1 : (m+m-k)*p-p );
end

