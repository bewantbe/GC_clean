% Give spike train from spike time data
% cal_mode == 0 means the spike train contains only 0 or 1 (time are rounded)
% cal_mode == 1 means the spike train is a average
% Can not handle multiple spikes at one time.

function [st, fr_cnt] = SpikeTrain(ras, len, neuron_id, bound_ignore, stv, cal_mode)
if ~exist('cal_mode', 'var')
  cal_mode = 0;
end
if ~exist('bound_ignore', 'var') || isempty(bound_ignore)
  bound_ignore = 1;
end
if ~exist('stv', 'var') || isempty(stv)
  stv = 0.5;
end

switch cal_mode
case 0
  acce_index = round(ras(ras(:,1)==neuron_id, 2) / stv);
  acce_index(acce_index >= len-bound_ignore | acce_index <= bound_ignore) = [];
  fr_cnt = length(acce_index);
  st = zeros(1,len);
  st(acce_index') = 1;
case 1  % 40% slower
  fire_id = ras(ras(:,1)==neuron_id, 2) / stv;   % exact fire time id
  acce_index = floor(fire_id);
  id_ignore = acce_index >= len-bound_ignore | acce_index <= bound_ignore;
  acce_index(id_ignore) = [];  fire_id(id_ignore) = [];
  fr_cnt = length(acce_index);
  fire_id = fire_id - acce_index;
  st = zeros(1,len);
  st(acce_index') = 1-fire_id;
  st(1+acce_index') = fire_id;
otherwise
  error('no this case');
end

end
