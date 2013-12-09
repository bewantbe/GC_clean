% Move the frequency 0 to the center
% i.e. from [0 1 2 3] to [-2 -1 0 1]
% S is len*p*p matrix

[S_np, fq_np] = fftFreq2seqFreq(S, freq)

N     = length(freq);
df    = freq(2)-freq(1);
fq_np = [ -(ceil((N-1)/2):-1:1)*df 0 (1:floor((N-1)/2))*df ];
S_np  = fftshift(S,1);

