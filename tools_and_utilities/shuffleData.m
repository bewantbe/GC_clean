% shuffle data

function X = shuffleData(X)
[p,len] = size(X);

for k=1:p
  [~,dd] = sort(rand(1,len));
  X(k,:) = X(k,dd);
end
%X = X'
%X = X(randperm(len,p));
%X = X';
end
