%

noisecov = diag([0.3 1.0 0.2]);
A = [-0.8  0.0  0.0  0.5  0.0  0.0;
      0.2 -0.9  0.0 -0.1  0.8  0.0;
      0.1  0.0 -0.5 -0.2  0.0  0.2];

fftlen = 1024;
S = A2S(A, noisecov, fftlen);
m = 50;
R = S2cov(S, m);

gc = RGrangerT(R)

tic;
len=1e7;
p = size(A,1);
%X = gdata(len, m, p, noisecov, A);
save('-ascii','-double', 'a_coef.txt', 'A');
halfde = chol(noisecov)';
save('-ascii','-double', 'd.txt', 'halfde');
tic;
system(['./GCcal/gendata -half-noise d.txt -l ',num2str(len),' -seed ',num2str(2^31*rand())]);
toc;
%output: out.txt raw file
fid = fopen('out.txt', 'r');
X = fread(fid, [p, Inf], 'double');
fclose(fid);
toc;

%gc_x = pos_nGrangerT2(X, m)

%gc_x2 = pos_nGrangerT2([X(1,:); X(2,:)+X(3,:)], m)

gc_x2 = pos_nGrangerT2([X(1,:); X(2,:)], m)

x = X(1,:);
y = X(2,:);
gc_x2 = pos_nGrangerT2([x-x.^3+x.^5; y.^3], m)

