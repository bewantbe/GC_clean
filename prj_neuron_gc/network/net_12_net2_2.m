%#!/usr/local/bin/octave -qf

pp = 12;

nnet = [
0 0 ones(1, pp-2);
1 0 ones(1, pp-2);
zeros(pp-2, pp);
];

fd = fopen('net_12_net2_2.txt','w');
nnet = int2str(int32(nnet));
for jj=1:pp
  fprintf(fd, '%s\n', nnet(jj,:));
end
fclose(fd);

