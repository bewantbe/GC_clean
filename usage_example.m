% GC_clean usage example.

% Addpath and check compilation of mex file.
% You may put it in startup.m
run('~/matcode/GC_clean/startup.m'); 

% Generate X(t) based on linear model (vector AutoRegression)
%   X(t) + A1*X(t-1) + A2*X(t-2) +...+ Am*X(t-m) = noise(t)

% Variance of noise(t)
De = diag([0.3 1.0 0.2]);

% Coefficients (always in concatenated form) A = [A1, A2]
A = [-0.8  0.0 -0.4  0.5 -0.2  0.0;
      0.0 -0.9  0.0  0.0  0.8  0.0;
      0.0 -0.5 -0.5  0.0  0.0  0.2];

len = 1e5;
X = gendata_linear(A, De, len);

% Choose a fitting order for GC.
od_max = 100;
[od_joint, od_vec] = chooseOrderFull(X, 'AICc', od_max);
m = max([od_joint, od_vec]);
% For a fast schematic test, use:
% m = chooseOrderAuto(X, 'AICc')

% The Granger Causality value (in matrix form).
GC = nGrangerTfast(X, m)

% Significance test: Non-zero probably based on 0-hypothesis (GC==0).
p_nonzero = gc_prob_nonzero(GC, m, len)

% The connectivity matrix.
p_value = 0.0001;
net_adjacency = p_nonzero > 1 - p_value

% This should give the same result.
gc_zero_line = chi2inv(1-p_value, m)/len;
net_adjacency2 = GC > gc_zero_line;

% Frequency domain GC.
fftlen = 1024;
wGC = nGrangerF(X, m, fftlen);
w = (0:fftlen-1)/fftlen;

figure(1);
p = size(X, 1);
for ii = 1 : p
  for jj = 1 : p
    if ii == jj continue; end
    subplot(p, p, (ii-1)*p + jj);
    plot(w, squeeze(wGC(ii, jj, :)), 'b', ...
         w, gc_zero_line*ones(size(w)), 'r');
  end
end

