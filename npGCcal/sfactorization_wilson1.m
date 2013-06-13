%any dim
% Using wilson's algorithm, we factorize spectral matrix to obtain transfer
% function and error covariance matrix.
% INPUT: S, spectral density matrix for 2 channels(channels x channels x
% frequency),
% freq, an array of frequency at which spectral densities are evaluated.
% fs, sampling frequency.
% OUTPUT: Snew, improved spectral density matrix (channels x channels x frequency),
% Hnew: transfer function(channels x channels x frequency),
% Znew: Error covariance matrix(channels x channels ),
function [Snew,Hnew,Znew] = sfactorization_wilson1(S,freq)
%Wilsons algorithm: 
m = size(S,1);N=length(freq)-1; tol = 1E-5;
f_ind=0;
for f=freq
    f_ind=f_ind+1;
    Sarr(:,:,f_ind)=S(:,:,f_ind);
    if(f_ind>1)
        Sarr(:,:,2*N+2-f_ind)=S(:,:,f_ind).';%'.0'
    end
end
% perform ifft to obtain gammas

for k1=1:m
    for k2=1:m
        gam(k1,k2,:)=real(ifft(squeeze(Sarr(k1,k2,:))));
    end
end

gam0 = gam(:,:,1);
h = chol(gam0); % this is psi_1
%h = gam0;
for ind = 1: size(Sarr,3),
       psi(:,:,ind) = h; % initialization for the 1st iteration
end

I = eye(m); Niterations = 100;

for iter = 1: Niterations,
       NOTEITE=iter;
      for ind = 1: size(Sarr,3),
            g(:,:,ind)=inv(psi(:,:,ind))*Sarr(:,:,ind)*inv(psi(:,:,ind))'+I;%'
     end

     gp = PlusOperator(g,m,freq);
psiold=psi;
   for k = 1: size(Sarr,3),
          psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
          psierr(k)=norm(psi(:,:,k)-psiold(:,:,k),1);
   end
   iter;
   psierrf=mean(psierr);
   if(psierrf<tol),break;end;
     
     
     %set error check & break statement here
   
end 

NOTEITE;
for k = 1: length(freq),
      Snew(:,:,k) = squeeze(psi(:,:,k)*psi(:,:,k)'); %'
end

for k1=1:m
    for k2=1:m
        gamtmp(k1,k2,:)=real(ifft(squeeze(psi(k1,k2,:))));
    end
end
A0=squeeze(gamtmp(:,:,1)); % this is psi_1
A0inv=inv(A0);
Znew=A0*A0.';
for k = 1: length(freq),
      Hnew(:,:,k) = squeeze(psi(:,:,k))*A0inv; %
      Serr(:,:,k)=Sarr(:,:,k)-Hnew(:,:,k)*Znew*Hnew(:,:,k)';
      Serrnorm(k)=norm(Serr(:,:,k),1)/norm(Sarr(:,:,k));
end


% figure(1)
% subplot(2,1,1)
% plot(freq,squeeze(Snew(1,1,:)),'g',freq,squeeze(S(1,1,:)),'b'); legend('Snew','S');
% subplot(2,1,2)
% plot(freq,Serrnorm); legend('Difference');
%---------------------------------------------------------------------
function gp = PlusOperator(g,m,freq);

%This function is for [ ]+operation

for k1=1:m
    for k2=1:m
          gam(k1,k2,:)= ifft(squeeze(g(k1,k2,:)));
    end
end

% take positive lags only and half of the zero lag

gamp = gam;beta0 = 0.5*gam(:,:,1); 
gamp(:,:,1) = triu(beta0);  %this is Stau
gamp(:,:,length(freq)+1:end) = 0;

% reconstitute

for k1=1:m
    for k2=1:m
         gp(k1,k2,:)= fft(squeeze(gamp(k1,k2,:)));
    end
end