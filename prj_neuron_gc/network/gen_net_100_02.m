% generate big 0/1 random netowrk.

netstr = 'net_100_02';

rand('state',5);            % set initial state
prob_one = 0.2;
p=100;
nnet = rand(p,p)<prob_one;
nnet(eye(p)==1) = 0;

nnet = nnet .* rand(p);
fd = fopen([netstr,'.txt'],'w');
nnet = num2str(nnet);
for jj=1:p
  fprintf(fd, '%s\n', nnet(jj,:));
end
fclose(fd);

