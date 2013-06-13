% Here, we evaluate Granger causality from spectral density matrix,
% transfer function and the error covariance matrix.
% INPUT: Spectral denstiy matrix for 2 channels Snew (channels x channels x frequency),
% Transfer function Hnew (channels x channels x frequency),
% Error covariance matrix Znew (channels x channels),
% frequency f at which causalities are evaluated,
% Sampling frequency fs.

% OUTPUT: Fx2y : Granger causality from 1st channel to second channel,
% Fy2x: Granger causality from 2nd channel to 1st channel,
% pp1: power spectrum of 1st channel,
% pp2: power spectrum of 2nd channel.
function [Fx2y,Fy2x,pp1,pp2]=pwcausalr_c(Snew,Hnew,Znew,freq)
L=2;
i=1;
j=2;
index=1;
f_ind = 0;
         for f = freq
             f_ind = f_ind + 1;
             pp1(f_ind) = abs(Snew(1,1,f_ind)*2);
             %if (i==L-1)&(j==L)
                 pp2(f_ind) = abs(Snew(2,2,f_ind)*2);
            %end
             eyx = Znew(2,2) - Znew(1,2)^2/Znew(1,1); %corrected covariance
             exy = Znew(1,1) - Znew(2,1)^2/Znew(2,2);
              Iy2x(index,f_ind) = abs(Hnew(1,2,f_ind))^2*eyx/abs(Snew(1,1,f_ind)); %measure within [0,1]
              Ix2y(index,f_ind) = abs(Hnew(2,1,f_ind))^2*exy/abs(Snew(2,2,f_ind));
             Fy2x(index,f_ind) = log(abs(Snew(1,1,f_ind))/abs(Snew(1,1,f_ind)-(Hnew(1,2,f_ind)*eyx*conj(Hnew(1,2,f_ind))))); %Geweke's original measure
             Fx2y(index,f_ind) = log(abs(Snew(2,2,f_ind))/abs(Snew(2,2,f_ind)-(Hnew(2,1,f_ind)*exy*conj(Hnew(2,1,f_ind)))));
             %Fxy(index,f_ind) = log(abs(Snew(1,1,f_ind)-(Hnew(1,2,f_ind)*eyx*conj(Hnew(1,2,f_ind)))/fs)*abs(Snew(2,2,f_ind)-(Hnew(2,1,f_ind)*exy*conj(Hnew(2,1,f_ind)))/fs)/abs(det(Snew(:,:,f_ind)));
         end
