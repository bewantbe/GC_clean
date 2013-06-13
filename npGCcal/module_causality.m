%any dim
% We have the spectrum evaluated from the binned data using multitaper 
% estimation or any other method. Using these estimates, we will find transfer function
% and error covariance matrix using Wilson's algorithm. Then evaluate
% Granger causality for different pairs of channels.
%INPUT:- Sp: spectral density matrix where Sp(i,j,k) is cross spectral
%density between ith and jth channel/neuron at kth frequency.
% f: array containing frequency at which spectrums are evaluated.
% dt: Sampling period in seconds.
% OUTPUT:-  causality is a matrix such that causality(i,j)contains Granger causality
% from ith channel/neuron to jth channel/neuron summed over all frequency. Note that for
% convenience we put -1 at all diagonal entries.

function [causality]= module_causality(Sp,f)
[s1,s2,s3]=size(Sp);
% Sampling frequency
%Fs = 1/dt;
ft = 1;
In = 0; %for image num
for i=1:(s1-1)
    causality(i,i)=0;
    for j=i+1:s1
        [Snew,Hnew,Znew] = sfactorization_wilson1(Sp([i,j],[i,j],:),f(1:s3));
        [Fx2y,Fy2x,pp1,pp2]=pwcausalr_c(Snew,Hnew,Znew,f(1:s3));
        causality(i,j)= mean(Fx2y); 
        causality(j,i)= mean(Fy2x);
        CauFre(ft,:)  = Fx2y;
        CauFre(ft+1,:)= Fy2x;        
        ft=ft+2;
%         figure(2);
%         In = In+1;
%         subplot(s1*(s1-1)/2,2,In);
%         plot(f(1:s3),Fx2y);
%         title([num2str(i) ' to ' num2str(j)]);
%         In = In+1;
%         subplot(s1*(s1-1)/2,2,In);
%         plot(f(1:s3),Fy2x);
%         title([num2str(j) ' to ' num2str(i)]);
    end
end
causality(s1,s1)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot fre cau
% w = 1:256;
% figure(3);
% subplot(3,2,1);
% plot(w,CauFre(1,:));
% title('x->y cau');
% subplot(3,2,2);
% plot(w,CauFre(2,:));
% title('y->x cau');
% subplot(3,2,3);
% plot(w,CauFre(5,:));
% title('y->z cau');
% subplot(3,2,6);
% plot(w,CauFre(4,:));
% title('z->y cau');
% subplot(3,2,5);
% plot(w,CauFre(3,:));
% title('x->z cau');
% subplot(3,2,4);
% plot(w,CauFre(6,:));
% title('z->x cau');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%