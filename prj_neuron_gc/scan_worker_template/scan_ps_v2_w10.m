% scan: fix ps

% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat

% to distinguish different parallel program instances (also dir)
%signature = 'data_scan_ps/w_01_st';
signature = 'data_scan_ps/v2_w10';

% scan value sets
s_net  = {'net_2_2'};
s_time = [1e7];
s_scee = [0.01];
s_prps = 0.01*logspace(log10(0.32),log10(0.5), 30);
s_ps   = 0.01;  % pr=0.31:ISI=16000, pr=0.34:ISI=5200,  pr=0.5:ISI=180
s_stv  = 1/2;

maxod  = 99;
s_od   = 1:maxod;
hist_div = 0:20:10000;          % ISI

extst = '--RC-filter -q';

scan_ps_v2_loop;
