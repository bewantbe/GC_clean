% Check the roots of AR process (confirm that AR is stable)

function la = ARroots(A, b_show_plot)

[p, mp] = size(A);
B = zeros(mp, mp);
B(1:p, :) = -A;
B(p+1:mp+1:mp*(mp-p)) = 1;
la = eig(B);
if exist('b_show_plot','var') && b_show_plot
  if isreal(la)
    la = complex(la);
  end
  %compass(la);
  plot(la,'*');
  hold on
  plot(cos(0:0.01:2.01*pi), sin(0:0.01:2.01*pi), 'r');
  axis([-1.1,1.1,-1.1,1.1], 'square');
  hold off
end
r = max(abs(la));
if r >= 1
    warning('GC:AR:unstable', 'unstable process !!');
end

end
