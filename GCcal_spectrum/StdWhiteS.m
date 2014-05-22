% Whiten the auto spectrum of each variable

function WS = StdWhiteS(S)
if size(S,2)~=size(S,3)
  error('S shoule be fftlen*p*p matrix');
end

p   = size(S,2);
len = size(S,1);
X = zeros(len, p);
for k=1:p
  X(:,k) = S2X1D(real(S(:,k,k)));
end

WS = zeros(size(S));
for k1=1:p
  % later, we may ignore WS(:, k, k) part
  for k2=1:p
    WS(:, k1, k2) = (1./X(:,k1)) .* S(:, k1, k2) .* conj(1./X(:,k2));
  end
end

% ensure absolute real number
for k=1:p
  WS(:, k, k) = real(WS(:, k, k));
end

%WS(:,1,1) = real((1./X) .* S(:,1,1) .* conj(1./X));
%WS(:,1,2) = (1./X) .* S(:,1,2) .* conj(1./Y);
%WS(:,2,1) = (1./Y) .* S(:,2,1) .* conj(1./X);
%WS(:,2,2) = real((1./Y) .* S(:,2,2) .* conj(1./Y));
