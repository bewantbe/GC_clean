% generate big 0/1 random netowrk. 20% edges

netstr = 'net_100_01';

rand('state',5);            % set initial state
prob_one = 0.2;
p=100;
nnet = rand(p,p)<prob_one;
nnet(eye(p)==1) = 0;

fd = fopen([netstr,'.txt'],'w');
nnet = int2str(int32(nnet));
for jj=1:p
  fprintf(fd, '%s\n', nnet(jj,:));
end
fclose(fd);

