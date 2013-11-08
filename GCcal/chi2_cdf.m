% Cumulative distribution function of noncentral chi-squared distribution

function p = chi2_cdf(x,k,c)
  x(x<0) = 0;
  c(c<0) = 0;

  if ~exist('c','var') || c==0
    p = chi2cdf(x,k);
    return
  end

% A good approximation
  h = 1 - 2/3 * (k+c)*(k+3*c)/(k+2*c)^2;
  p = (k+2*c)/(k+c)^2;
  m = (h-1)*(1-3*h);
  v = ((x/(k+c)).^h - (1+h*p*(h-1-0.5*(2-h)*m*p))) / (h*sqrt(2*p)*(1+0.5*m*p));
  p = normcdf(v);

% By def
%  p = quadgk(@(x)chi2_pdf(x,k,c), 0, x);
% or
%  p = quad(@(x)chi2_pdf(x,k,c), 0, x);

% Exact but inf and might unstable
%{
  v = zeros(size(x));
  for jj = 0:100
    tmp_c = (c/2)^jj ./ gamma(jj+1);
    if (max(tmp_c) < 1e-16)
      break;
    end
    v = v + tmp_c .* chi2cdf(x, k+2*jj);
  end
  p = exp(-c/2) * v;
%}

end
