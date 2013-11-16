% generate big regular netowrk

netstr = 'net_100_10';

p=100;
n=p;
adj = zeros(p);
for k=1:2
  adj = adj + diag(ones(1,n-k),k);
  adj = adj + diag(ones(1,k),n-k);
  adj = adj + diag(ones(1,n-k),-k);
  adj = adj + diag(ones(1,k),k-n);
end
nnet=adj;

fd = fopen([netstr,'.txt'],'w');
nnet = num2str(nnet);
for jj=1:p
  fprintf(fd, '%s\n', nnet(jj,:));
end
fclose(fd);

