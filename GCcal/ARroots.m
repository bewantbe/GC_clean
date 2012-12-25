% Check the roots of AR process (confirm that AR is stable)

function la = ARroots(A)

[p, mp] = size(A);
B = zeros(mp, mp);
B(1:p, :) = -A;
B(p+1:mp+1:mp*(mp-p)) = 1;
la = eig(B);
la = la + 1e-300*I;            % convert all numbers to complex
%compass(la);
plot(la,'*', cos(0:0.01:2.01*pi), sin(0:0.01:2.01*pi));
axis([-1.1,1.1,-1.1,1.1], 'square');
r = max(abs(la));
if r >= 1
    disp('Warning: unstable process !!');
end

end
