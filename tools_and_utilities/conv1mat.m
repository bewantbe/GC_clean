% matrix conv
% only positive part
% A, B and C are in the 2-D format

function C = conv1mat(A, B)
if (size(A,1) ~= size(B,1))
    error('size(A,1) != size(B,1)');
end

p = size(A,1);
m = round(size(A,2)/p);
AA = zeros(p,p,(m+1)*p);
AA(:,:,1) = eye(p);
for k=1:m
    AA(:,:,k+1) = A(:, 1+(k-1)*p:k*p);
end
A = AA;

m = round(size(B,2)/p);
BB = zeros(p,p,(m+1)*p);
BB(:,:,1) = eye(p);
for k=1:m
    BB(:,:,k+1) = B(:, 1+(k-1)*p:k*p);
end
B = BB;

p  = size(A,1);
la = size(A,3);
lb = size(B,3);
l  = la+lb-1;
C  = zeros(p,p,l);
for tau=1:l
    c = zeros(p);
    for t=1:la
        t2 = tau-t+1;
        if ((t2>lb) || (t2<1)) % TODO: slow here
            continue;
        end
        c = c + A(:,:,t)*B(:,:,t2);
    end
    C(:,:,tau) = c;
end

CC = zeros(p,l*p-p);
CC(:,1:2) = C(:,:,1);
for k=1:l-1
    CC(:,1+(k-1)*p:k*p) = C(:,:,k+1);
end
C = CC;

end
