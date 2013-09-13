% test accuracy

%len = 2000;
%p = 2;
%n_trials = 30;
%x1 = randn(len, p, n_trials);
%t1 = rand(len, p)*2*pi;

%fq_max = 2000/2;
%fqs = (1:round(2*fq_max)) - floor(round(2*fq_max)/2) - 1;  % desired frequency points
%fqs = ifftshift(fqs);

%tic();
%xfa=nufftw(x1(:,1,2), t1(:,2), fq_max);
%toc()

%tic();
%err = zeros(1,p);
%for k=1:p
  %xf  = nuftw(x1(k,:), t1(k,:), fqs);
  %err(k) = 2*norm(xf - xfa(k,:))/sqrt(len);
  %err(k)
%end
%toc();

%mean(err)

%figure(3);
%hd=plot(abs(xf - xfa(end,:)), '-o');  set(hd, 'markersize', 2);

len = 2000;
p = 2;
n_trials = 30;
mX = randn(p, len, n_trials);
mT = rand(len, n_trials)*2*pi;
fftlen = len;
aveS = mX2S_nuft(mX, mT, fftlen);

