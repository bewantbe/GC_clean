%plot in arbitrary x scale
% fxs should be a monotonic function
% inv_fxs is inverse function of fxs
% if label_num is a vector, interpret it as the labels(ticks).
% if y is empty ([]), do not draw the plot, but return the x position, ticks position, tick values

function [H, x_tick, x_tick_val] = xscale_plot(x,y,fxs,inv_fxs, label_num)
  xl = fxs(x);
  xl_bg = xl(1);
  xl_ed = xl(end);
  if ~exist('label_num','var')
    label_num = 10;
  end
  if length(label_num) == 1
    x_tick = linspace(xl_bg, xl_ed, label_num);
  else
    x_tick = fxs(label_num);
  end
  x_tick_val = inv_fxs(x_tick);
  % shift to an easy to read number
  o = 10.^(floor(log10(x_tick_val))-1);
  %x_tick_val = round(x_tick_val./o).*o;
  % take special care of boundary?
  x_tick_val(2:end-1) = round(x_tick_val(2:end-1)./o(2:end-1)).*o(2:end-1);
  if xl_bg < xl_ed
    x_tick_val(1)   = ceil (x_tick_val(1)  ./o(1)  ).*o(1);
    x_tick_val(end) = floor(x_tick_val(end)./o(end)).*o(end);
  else
    x_tick_val(1)   = floor(x_tick_val(1)  ./o(1)  ).*o(1);
    x_tick_val(end) = ceil (x_tick_val(end)./o(end)).*o(end);
  end
  x_tick_val(o==0) = 0;
  x_tick = fxs(x_tick_val);
  if ~isempty(y)
    H = plot(xl, y);
    % draw x ticks (scale label)
    set(gca,'xtick', x_tick);
    set(gca,'xticklabel', x_tick_val);
  else
    H = xl;
  end
end
