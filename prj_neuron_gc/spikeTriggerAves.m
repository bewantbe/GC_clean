% spike trigger average

%r   : the guy who is firing
%s   : used to do average (allow to be array)
%ras : firing time
%X   : used to do averaging

% usage example:
% [X, ISI, ras] = gendata_neu('net_2_2', 0.02, 1, 0.012, 1e6, 0.5);
% [tg_ave_volt, s_rel_time] = spikeTriggerAve(1,2,ras,X, 200, 0.5);
% plot(s_rel_time, tg_ave_volt, '-+');

function [tg_ave, s_rel_time] = spikeTriggerAve(r, s, ras, X, ana_len, stv)
if (exist('ana_len','var')==0)
    ana_len = 200;         % length of analysis window
end
if (exist('stv','var')==0)
    stv = 0.5;             % sampling interval
end

if length(ana_len)==1
  % relative time offset range, in sample point unit
  s_rel_time = ceil(-ana_len/4+1) : ceil(ana_len/4*3);
else
  s_rel_time = ana_len(1):ana_len(2);
end

fire_index = round(ras(ras(:,1)==r, 2) / stv);
% delete head and tail for safe
len = size(X,2);
fire_index(fire_index > len-s_rel_time(end) | fire_index <= -s_rel_time(1)) = [];

tg_ave = zeros(length(s), length(s_rel_time));
for id_t0 = 1:length(s_rel_time)
    tg_ave(:, id_t0) = mean(X(s, fire_index+s_rel_time(id_t0)), 2);
end

s_rel_time = s_rel_time*stv;  % convert to time unit
end
