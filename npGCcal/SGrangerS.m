function [cdcausal, Znew, Hnew] = SGrangerS(Sp)
%format long;
nfft=2^10;
if (exist('Sp','var')==0)
	[Sp,R]=para2cov(2,3,nfft);
end
dt =0.01;
fs = 1/dt;
[s1,s2,s3]=size(Sp);
df=fs/s3;
f=0:df:fs;
% Evaluation of G. causality
M1=s1;
if M1 <3
%	[causality]= module_causality(Sp,f,fs);
	[causality]= module_causality(Sp,f);
	%disp('causality')
	causality=causality';
	%disp(causality)
	cdcausal = causality;
end
if M1 >=3
freq =f(1:s3);
% Evaluation of conditional G. causality
%tic();
%[Snew,Hnew,Znew] = sfactorization_wilson1(Sp,freq,fs);
[Snew,Hnew,Znew] = sfactorization_wilson1(Sp,freq);

i=0;
j=0;
k=0;
t=0;
 cdcausal=zeros(M1,M1);
while j <M1
    
    j=j+1;
    i = 0;
    while i <M1
        i = i+1;
        if i == j
            continue;
        end
        clear Sp2;
%         if j==1
%             Sp2=Sp(2:M1,2:M1,:);
%         end
%         if j>1 && j<M1
%             Sp2(1:j-1,1:j-1,:)=Sp(1:j-1,1:j-1,:);
%             Sp2(j:M1-1,j:M1-1,:)=Sp(j+1:M1,j+1:M1,:);         
%         end
%         if j==M1
%             Sp2=Sp(1:M1-1,1:M1-1,:);
%         end
        for ii=1:j-1
            for jj=1:j-1
                Sp2(ii,jj,:)=Sp(ii,jj,:);
            end
            for jj=j+1:M1
                Sp2(ii,jj-1,:)=Sp(ii,jj,:);
            end
        end
        for ii=j+1:M1
            for jj=1:j-1
                Sp2(ii-1,jj,:)=Sp(ii,jj,:);
            end
            for jj=j+1:M1
                Sp2(ii-1,jj-1,:)=Sp(ii,jj,:);
            end
        end
        if j>i
            t=i;
        end
        if j<i
            t=i-1;
        end
            [Snew2,Hnew2,Znew2] = sfactorization_wilson1(Sp2,freq);
%            [Snew2,Hnew2,Znew2] = sfactorization_wilson1(Sp2,freq,fs);
            cdcausal(i,j)= log(Znew2(t,t)/Znew(i,i));     
        end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:M1
    cdcausal(i,i)= 0;
end
%disp('Condn G.causality')
%% disp('the kth matrix mean the causality on condition kth chanel');
%disp(cdcausal)
%%save cdcausal;
%toc();
end
