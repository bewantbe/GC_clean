% Convert covariance series (function) to big covariance matrix
% like (R(i)=E( X(t)*X(t-i)' )):
% R = [R(0), R(1), R(2)];
% covz = R2covz(R);
% covz = 
%   [ R(0)  R(1)  R(2)
%     R(1)' R(0)  R(1)
%     R(2)' R(1)' R(0) ]

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

