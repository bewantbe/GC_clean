% Expample
%   [cf3d, cb3d, df3d, db3d, A] = BlockLevinson_v3_inv(R, 1);
%   QyyArr = GetDiagInvCov(cf3d, cb3d, df3d, db3d)

function QyyArr = GetDiagInvCov(cf3d, cb3d, df3d, db3d)
p = size(cf3d,1);
m = size(cf3d,3);

cf3d = permute(cf3d, [2,1,3]);  % use the transposed form
df3d = permute(df3d, [2,1,3]);  % use the transposed form

QyyArr = zeros(p,m,m);
f_T2 = fft(db3d(:,:,[end, 1:end-1]), m, 3);
f_T4 = fft(cb3d(:,:,[end, 1:end-1]), m, 3);
%f_T2 = fft(cat(3,db3d(:,:,[end, 1:end-1]),zeros(p,p,m==1)), m, 3);
%f_T4 = fft(cat(3,cb3d(:,:,[end, 1:end-1]),zeros(p,p,m==1)), m, 3);
df3d(:,:,1) = df3d(:,:,1) - eye(p);
% The FFT in Octave(3.6.2) can not do fft along singleton 3rd dimension.
% Like fft(rand(2,2,1), 9, 3)
j = 1;
f_T1 = fft(cat(3,cf3d(:,:,1),zeros(p,p,1)), m, 3);
f_T3 = fft(cat(3,df3d(:,:,1),zeros(p,p,1)), m, 3);
QyyArr(:,:,j) = squeeze(...
  sum(f_T1.*f_T2, 1) -...
  sum(f_T3.*f_T4, 1));
QyyArr(:,:,j) = real(ifft(QyyArr(:,:,j), m, 2));

for j=2:m
  QyyArr(:,:,j) = real(ifft(squeeze(...
    sum(fft(cf3d(:,:,j:-1:1), m, 3).*f_T2) -...
    sum(fft(df3d(:,:,j:-1:1), m, 3).*f_T4)...
    ), m, 2));
end
QyyArr = permute(QyyArr, [3,2,1]);

