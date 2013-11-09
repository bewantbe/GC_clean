% scan: fix ps

% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat

% to distinguish different parallel program instances (also dir)
%signature = 'data_scan_ps/w_01_st';
signature = 'data_scan_ps/v2_w1';

% scan value sets
s_net  = {'net_2_2'};
s_time = [1e5];
s_scee = [0.02];
s_prps = 0.9*[0.005:0.00025:0.00725, 0.0075:0.0005:0.0095, 0.01:0.001:0.019, 0.02:0.002:0.04];
s_ps   = [0.01, 0.02];
%s_ps   = [0.005, 0.01, 0.02, 0.03];
%s_prps = 0.9*[0.005:0.0005:0.00725, 0.0075:0.001:0.0095, 0.01:0.001:0.019, 0.02:0.002:0.04];
%s_ps   = [0.01, 0.02];
%s_prps = [0.005, 0.0075, 0.01, 0.02];
%s_ps   = [0.01, 0.02];
s_stv  = [16]/32;

maxod  = 99;
s_od   = 1:maxod;
hist_div = 0:1:500;          % ISI

extst = '--RC-filter -q';

scan_ps_v2_loop;
