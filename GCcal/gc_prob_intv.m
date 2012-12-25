% Calculate confidence interval of GC.
% only suitable for nGrangerT series functions(use same fitting order
%  in auto and joint regression).

function [gc_lower, gc_upper] = gc_prob_intv(gc, od, len, a)
if nargin<3
  disp(' Calculate confidence interval of GC.');
  usage('[gc_lower, gc_upper] = gc_prob_intv(gc, od, len [, a])');
end
if (exist('a','var')==0)
  a = 0.95;
end

%gc = exp(gc)-ones(size(gc));     % for more accuracy, but actually seems worse

% 95% confidence interval, good for od>16, suggest by Geweke 1982
%za = sqrt(2)*erfinv(a*0.947);   % *0.947 is a correction obtained from experiment.
za = sqrt(2)*erfinv(a*0.97);
gc_lower = (sqrt(gc - (od-1)./(3*len)) - za./sqrt(len)).^2 - (2*od+1)./(3*len);
gc_upper = (sqrt(gc - (od-1)./(3*len)) + za./sqrt(len)).^2 - (2*od+1)./(3*len);
if (size(gc,1)>1 && size(gc,1)==size(gc,2))
  gc_lower(eye(size(gc))==1) = 0;
  gc_upper(eye(size(gc))==1) = 0;
end
% avoid negative or non-usual number
gc_lower(imag(gc_lower)~=0) = 0;
gc_lower(gc_lower<0) = 0;
gc_upper(imag(gc_upper)~=0) = 0;
gc_upper(gc_upper<0) = 0;

%gc_lower = log(ones(size(gc))+gc_lower);
%gc_upper = log(ones(size(gc))+gc_upper);

end
