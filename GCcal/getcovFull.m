% get time symmetric cov series
% usage:
%   cov_full = getcovFull(X, m);
% or
%   cov_full = getcovFull(R);

function [cov_full, idt] = getcovFull(X, m)
if ~exist('m', 'var')
  % only covariance is porvided
  R = X;
else
  R = getcovpd(X, m);
end
p = size(X,1);
m3R = reshape(R,p,p,[]);
cov_full = cat(3,flipdim(permute(m3R(:,:,2:end), [2 1 3]),3),m3R);
idt = 1-size(m3R,3) : size(m3R,3)-1;
