% generate group network from simple network

neu_network = load('-ascii', 'net_3_06.txt');
kp = 5;                     % neurons in each group
rand('state',10);            % set initial state
in_group_strength = 0.6;
out_group_strength = 0.3;

p=size(neu_network,1);
pp = p*kp;                  % total number of nerons
nnet = zeros(pp);
for j1=1:p
for j2=1:p
  if j1==j2
    tmp = (rand(kp)<in_group_strength) & (eye(kp)==0);
  else
    tmp = (rand(kp)<out_group_strength) * neu_network(j1,j2);
  end
  nnet((j1-1)*kp+1:j1*kp, (j2-1)*kp+1:j2*kp) = tmp;
end
end

fd = fopen('net_3g5_06.txt','w');
nnet = int2str(int32(nnet));
for jj=1:pp
  fprintf(fd, '%s\n', nnet(jj,:));
end
fclose(fd);

