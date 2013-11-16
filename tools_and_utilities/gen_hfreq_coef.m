% generate a high frequency AR coefficient

% G is coefficient without the leading "one"
% de is corresponding variance so that variance
%    of x is 1.0

% r is the absolute value of the roots
% f is the frequency, rangre: -0.5~0.5

% e.g.
%  [G, de] = gen_hfreq_coef(0.9, 0.05);
%  plot(real(squeeze(A2S(G, de, 1024)))); legend(['Order ',num2str(length(G))]);

function [G, de] = gen_hfreq_coef(r, f, od)
if ~exist('od','var')
    od = 8;
end
pl = [1, -r*exp(-i*2*f*pi)];
G2 = real(conv(pl,conj(pl)));
G = G2;
for k = 1:round(od/2-1)
    G = conv(G2,G);
end
G = G(2:end);
S = real(squeeze(A2S(G, 1, 1024)));
de = 1/mean(S);    % so that mean(S)==1

end
