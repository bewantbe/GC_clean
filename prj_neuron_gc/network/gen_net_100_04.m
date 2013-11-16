% generate big random netowrk, entries are normal distributed.

netstr = 'net_100_04';

rand('state',5);            % set initial state
prob_one = 0.2;
p=100;
nnet = rand(p,p)<prob_one;
nnet(eye(p)==1) = 0;

nnet = nnet .* (1+randn(p)/5);
nnet(nnet<=0) = 0;
fd = fopen([netstr,'.txt'],'w');
nnet = num2str(nnet);
for jj=1:p
  fprintf(fd, '%s\n', nnet(jj,:));
end
fclose(fd);

