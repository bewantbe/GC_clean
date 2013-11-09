
function [nn, xx] = hist_pdf(x, n_bin);

[nn,xx] = hist(x, n_bin);
nn = nn/length(x)/(xx(end)-xx(1))*length(nn);

