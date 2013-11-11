
len = 1e5;

p = 5;
[q,r] = qr(randn(p));
v = logspace(-8, 7, p);

M = q*diag(sqrt(v));
X = M*randn(p, len);

cond(X*X'/len)

od = 4;
tic();
[g,d] = pos_nGrangerT(X, od);
toc();
tic();
[g2,d2] = pos_nGrangerT2(X, od);
toc(); tic();
[g3,d3] = pos_nGrangerT_qr(X, od);
toc();
