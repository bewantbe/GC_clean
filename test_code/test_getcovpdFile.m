p=2; fn='data/volt_net_2_2_sc=0.01_pr=1_ps=0.012_t=1.01e+06_stv=0.5.dat';
p=100; fn='data/volt_net_100_21_p[80,20]_sc=[0.004,0.004,0.008,0.008]_pr=1_ps=0.012_t=1.00e+06_stv=0.5.dat';

tic; fid=fopen(fn,'r'); X=fread(fid,[p,inf],'double'); fclose(fid); toc
size(X)

m = 20;
tic; Ro = getcovpd(X, m); toc

tic; R = getcovpdFile(fn, p, m); toc

R - Ro

max(abs(R-Ro)(:))

%R = getcovpdFile(fn, 2, m, 20000);


%{

[p, len] = size(X);
aveX = mean(X, 2);
m = 2;
k = 0;

(X(:, m+1:end) - aveX) * (X(:,m+1-k:end) - aveX)' / (len-m)

Rt = X(:, m+1:end) * X(:,m+1-k:end)'/(len-m) - aveX*aveX'

Rt...
+ aveX*aveX'...
-aveX * sum(X(:,m+1-k:end),2)'/(len-m)...
+ aveX*aveX'...
-sum(X(:, m+1:end),2) * aveX'/(len-m)

Rt...
+ aveX*(aveX*len - sum(X(:,m+1-k:end),2) - m*aveX)'/(len-m)...
+ (aveX*len - sum(X(:, m+1:end),2) - m*aveX) * aveX'/(len-m)
%}
