% Solve block toeplitz matrix by Levinson algorithm
% Like other levinson-like algorithm, this program does not suitable for
% large condition number problem (say k<10^6)
% (this is BlockLevinson_v4_o2.m in the test_blocklevinson)

% Mathematically, we have
%   [Aall, Deps] = ARregression(R)
% Aall = af, s_de(:,:,end) = Deps

function [af, s_de] = BlockLevinson(R)
[p, m] = size(R);
m = round(m/p)-1;
if m<1
  error('no enough input data!');
end

b_R1 = reshape(flipdim(permute(reshape(...
         R(:, p+1:end),p,p,[]) ,[2,1,3]) ,3) ,p,[]);

af = - R(:,1:p) \ b_R1(:,end-p+1:end);
ab = - R(:,1:p) \ R(:,p+1:2*p);
df = R(:,1:p) + R(:,p+1:2*p) * af;
db = R(:,1:p) + b_R1(:,end-p+1:end) * ab;
s_de = zeros(p,p,m);
s_de(:,:,1) = df;

for k=2:m
  ef = b_R1(:, end-p*(k-1)+1:end) * af + b_R1(:, end-p*k+1:end-p*(k-1));
  eb = R(:, p+1:p*k) * ab + R(:, p*k+1:p*(k+1));
  gf = db \ ef;
  gb = df \ eb;
  tmp_af = af;
  af = [tmp_af; zeros(p)] - [ab; eye(p)] * gf;
  ab = [zeros(p); ab] - [eye(p); tmp_af] * gb;
  df = df - eb * gf;
  db = db - ef * gb;
  s_de(:,:,k) = df;
end

af = af';  % we need this form, in compatible with ARregression.m

end
