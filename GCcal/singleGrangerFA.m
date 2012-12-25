% Core of Frequency domain Granger causality
% input coefficients D, B and noise term De0_star, De0
% return wGc the GC from id_y to id_x with length fftlen

function wGc = singleGrangerFA(D, De0_star, B, De0, id_y, id_x, fftlen)
p = size(B,1);
lx = length(id_x);  ly = length(id_y);  lz = p-lx-ly;

% prepare variable indexes
id_z = 1:p;
id_z([id_x, id_y]) = [];    % rest veriables
tmprg = int16(zeros(1,p));
tmprg(id_x) = 1;
tmprg(id_z) = 3;
tmprg(id_y) = [];
id_d_x = find(tmprg==1);
id_d_z = find(tmprg==3);    % for autoregression

%% renormalization
% for autoregression
D0=eye(lx+lz);                                     % renormalizing matrix for x,z
D0(id_d_z, id_d_x) = -De0_star(id_d_z,id_d_x) / De0_star(id_d_x,id_d_x);
De_star = D0*De0_star*D0';                         % renormalized noise variance

% for joint-regression
Q1 = eye(p);
Q1(id_y,id_x) = -De0(id_y,id_x) / De0(id_x,id_x);
Q1(id_z,id_x) = -De0(id_z,id_x) / De0(id_x,id_x);
yDe = Q1*De0*Q1';
Q2 = eye(p);
Q2(id_z,id_y) = -yDe(id_z,id_y) / yDe(id_y,id_y);  % eliminate cov(y,z|x,y,z)
B0 = Q2*Q1;                                        % renormalizing matrix for x,y,z
De = B0*De0*B0';                                   % renormalized noise variance

%%GC
%sprintf('%19.16f', log(det(De0_star(id_d_x,id_d_x))/det(De0(id_x,id_x))))
%sprintf('%19.16f', log(det( De_star(id_d_x,id_d_x))/det( De(id_x,id_x))))

%% turn to frequency domain
% autoregressive coef
lD = cat(3, eye(size(D,1)), reshape(D, [size(D,1),size(D,1),size(D,2)/size(D,1)]));
for k=1:size(lD,3)
  lD(:,:,k) = D0*lD(:,:,k);
end
wD = fft(lD, fftlen, 3);

% moving average coef of joint-regression
lB = cat(3, eye(size(B,1)), reshape(B, [size(B,1),size(B,1),size(B,2)/size(B,1)]));
for k=1:size(lB,3)
  lB(:,:,k) = B0*lB(:,:,k);
end
wB = fft(lB, fftlen, 3);
wA = zeros(p,p,fftlen);
for k=1:fftlen
  wA(:,:,k) = inv(wB(:,:,k));
end

% at last...
wF11 = zeros(lx,lx,fftlen);
wF12 = zeros(lx,ly,fftlen);
wF13 = zeros(lx,lz,fftlen);
wGc  = zeros(1,fftlen);
% Use this (sxs) to get a exactly match to nGrangerT
% But might obtain a non-positive causality
%sxs  = det(De_star(id_d_x,id_d_x));
for k=1:fftlen
  wF11(:,:,k) = wD(id_d_x,id_d_x,k)*wA(id_x,id_x,k) + wD(id_d_x,id_d_z,k)*wA(id_z,id_x,k);
  wF12(:,:,k) = wD(id_d_x,id_d_x,k)*wA(id_x,id_y,k) + wD(id_d_x,id_d_z,k)*wA(id_z,id_y,k);
  wF13(:,:,k) = wD(id_d_x,id_d_x,k)*wA(id_x,id_z,k) + wD(id_d_x,id_d_z,k)*wA(id_z,id_z,k);
  sxs = det(wF11(:,:,k)*De(id_x,id_x)*wF11(:,:,k)' + wF12(:,:,k)*De(id_y,id_y)*wF12(:,:,k)' + wF13(:,:,k)*De(id_z,id_z)*wF13(:,:,k)');
  wGc(k) = log(real(sxs/det(wF11(:,:,k)*De(id_x,id_x)*wF11(:,:,k)')));
end

end
