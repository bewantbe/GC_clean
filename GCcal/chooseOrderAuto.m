% Choose Order by AIC/BIC
% if od_full not provided, stop when best order is found.
%   otherwise, calculate all BIC/AIC till the order "od_full".
% Try levinson algorithm first. if the condition number seems very large,
%   use positive-defined covz to obtain the order.

function [xic_od, s_xic_val, s_lndet_de] = chooseOrderAuto(X, ic_mode, od_full)
b_lack_accuracy = false;
if ~exist('od_full', 'var')
  od_max = 500;
else
  if od_full<0
    b_lack_accuracy = true;
  end
  od_max = abs(od_full);
end
if ~exist('ic_mode', 'var') || isempty(ic_mode)
  ic_mode = 'BIC';
end
logdet = @(B) 2*sum(log(diag(chol(B))));
[p,len] = size(X);
ic_mode = upper(ic_mode);
switch ic_mode
  case 'BIC'
    fxic = @(lde,m) lde + p^2*m*log(len)/len;
  case 'AIC'
    fxic = @(lde,m) lde + 2*p^2*m/len;
  case 'AICC'
    fxic = @(lde,m) lde + 2*p^2*m/len + 2*p^2*m*(p^2*m+1)/(len-p^2*m-1)/len;
  otherwise
    error('ic_mode: only accept BIC,AIC,AICc');
end

if ~b_lack_accuracy  % high accuracy is not requested, so try Levinson first
  %disp('Trying Levinson...');
  X = bsxfun(@minus, X, mean(X,2));
  Rxy0 = (X(:, od_max+1:len) * X(:, od_max-0+1:len-0)') / (len-od_max);
  Rxy1 = (X(:, od_max+1:len) * X(:, od_max-1+1:len-1)') / (len-od_max);
  R1 = Rxy1;
  b_R1 = Rxy1';

  af = - Rxy0 \ b_R1;
  ab = - Rxy0 \ Rxy1;
  df = Rxy0 + Rxy1 * af;
  db = Rxy0 + b_R1 * ab;
  s_lndet_de = logdet(df);
  xic_min = fxic(s_lndet_de(end),1);
  xic_od  = 1;
  mono_len = 0;
  xic_tmp_old  = xic_min;
  s_xic_val = xic_min;
  diag_Rxy0 = diag(Rxy0);

  for k = 2:od_max
    Rxy = (X(:, od_max+1:len) * X(:, od_max-k+1:len-k)') / (len-od_max);

    ef = b_R1 * af + Rxy';
    b_R1 = [Rxy', b_R1];
    eb = R1 * ab + Rxy;
    R1 = [R1, Rxy];
    gf = db \ ef;
    gb = df \ eb;
    tmp_af = af;
    af = [tmp_af; zeros(p)] - [ab; eye(p)] * gf;
    ab = [zeros(p); ab] - [eye(p); tmp_af] * gb;
    df = df - eb * gf;
    db = db - ef * gb;

    % in case of high condition number, try the slow but accurate method
    if any(diag(lu(df))<=0) || min(diag(df)./diag_Rxy0) < sqrt(0.5/len)
      b_lack_accuracy = true;
      %min(diag(df)./diag_Rxy0)
      %diag(df)./diag_Rxy0
      warning('chooseOrderAuto: lack of accuracy, trying SPD covz ...');
      break;
    end
    s_lndet_de = [s_lndet_de, logdet(df)]; 
    xic_tmp = fxic(s_lndet_de(end),k);
    s_xic_val = [s_xic_val, xic_tmp];
    % judge best order
    if xic_min > xic_tmp
      xic_min = xic_tmp;
      xic_od  = k;
    end
    if ~exist('od_full','var')
      % try some more fit, to make sure it's minimum (may fail in extreme case)
      if (xic_tmp_old < xic_tmp)
        mono_len = mono_len + 1;
        if mono_len >= 10
          %tried_until_od = k
          break;
        end
      else
        mono_len = 0;
      end
      xic_tmp_old = xic_tmp;
    end
  end  % for k
end

if ~b_lack_accuracy
  % solution obtained successfully
  return
end

%disp('Using covz version...');
s_xic_val = [];
s_lndet_de = [];
if exist('k', 'var') && ~exist('od_full','var')  % if we failed in Levinson
  od_max = min([2*k+10, od_max]);
  warning(sprintf('chooseOrderAuto: od_max set to %d. Specify the order explicitly if you want something different.', od_max));
end
covz = getcovzpd(X, od_max);
for k = 1:od_max
  [~, De] = ARregressionpd(covz(end-(k+1)*p+1:end, end-(k+1)*p+1:end), p);
  s_lndet_de = [s_lndet_de, logdet(De)]; 
  xic_tmp = fxic(s_lndet_de(end),k);
  s_xic_val = [s_xic_val, xic_tmp];
end
[xic_min, xic_od] = min(s_xic_val);

end
