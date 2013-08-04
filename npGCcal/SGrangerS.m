function [cdcausal, Znew, Hnew] = SGrangerS(Sp)
Sp = ipermute(Sp, [3,1,2]);   % convert to p*p*fftlen
[p,s2,mlen]=size(Sp);
f = 0:1/mlen:1-1/mlen;

% Evaluation of G. causality
if p <3
    [causality]= module_causality(Sp,f);
    cdcausal=causality';
end

if p >=3
freq =f(1:mlen);
% Evaluation of conditional G. causality
[Snew,Hnew,Znew] = sfactorization_wilson1(Sp,freq);

i=0;
j=0;
k=0;
t=0;
cdcausal=zeros(p,p);
while j <p
    
    j=j+1;
    i = 0;
    while i <p
        i = i+1;
        if i == j
            continue;
        end
        for ii=1:j-1
            for jj=1:j-1
                Sp2(ii,jj,:)=Sp(ii,jj,:);
            end
            for jj=j+1:p
                Sp2(ii,jj-1,:)=Sp(ii,jj,:);
            end
        end
        for ii=j+1:p
            for jj=1:j-1
                Sp2(ii-1,jj,:)=Sp(ii,jj,:);
            end
            for jj=j+1:p
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
%         [Snew2,Hnew2,Znew2] = sfactorization_wilson1(Sp2,freq,fs);
        cdcausal(i,j)= log(Znew2(t,t)/Znew(i,i));     
        end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:p
    cdcausal(j,j)= 0;
end
end
