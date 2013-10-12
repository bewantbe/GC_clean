%

% _ST _expIF
signature = 'data_scan_stv/IF_w53';

% scan value sets
s_net  = {'net_2_2'};
s_time = [5e7];
s_scee = [0.02];
s_prps = [0.012];
s_ps   = [0.012];
s_stv  = [0.5:0.5:20];
maxod  = 99;
s_od   = 1:maxod;
hist_div = 0:0.5:400;         % ISI
T_segment = 1000;             % in ms
stv0   = 0.125;               % fine sample rate

resample_mode = 'r';
b_overlap_time_interval = true;

scan_stv_work_50_loop;

