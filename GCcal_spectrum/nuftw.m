% naive nuft

function Xf = nuftw(X,T,fqs, md)

if ~exist('md','var') || md == 0

%Xf = zeros(size(fqs));
%for j=1:length(T)
  %Xf += exp(-I*T(j)*fqs)*X(:,j);
%end

Xf = zeros(size(fqs));
for k=1:length(fqs)
  Xf(k) = exp(-I*T*fqs(k))*X';
end

else

% use Kahan summation
Xf = zeros(size(fqs));
c  = zeros(size(fqs));
for j=1:length(T)
  y = exp(-I*T(j)*fqs)*X(:,j) - c;
  t = Xf + y;
  c = (t - Xf) - y;
  Xf = t;
end

end
