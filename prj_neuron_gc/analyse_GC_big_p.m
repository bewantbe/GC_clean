% analysis single GC trial with large p
% run analyse_GC_simple.m first to get data.
is_octave = exist('OCTAVE_VERSION','builtin') ~= 0;
if is_octave
%    font_size = 28;
    font_size = 24;
    line_width = 5;
else
    font_size = 24;
    line_width = 2;
end

pic_prefix0 = 'pic_tmp/';

set(0, 'defaultfigurevisible', 'off');

auto_gc_zero_cut = false;
b_output_pics    = true;
b_cal_net        = false;
b_use_pairGC     = false;

p_val = 2e-4;  % aic_od 1e-15 and 1.65 for S=0.01, bic_od 1e-2 and 1.015 for S=0.005

s_case_th = [117];  %[20 2 21 22 23 24] %[1 2 3 4 5 6 7 8]  [32 33 34]
n_case = length(s_case_th);

s_stv= zeros(1, n_case);
s_aic= zeros(1, n_case);
s_bic= zeros(1, n_case);
s_CC = zeros(1, n_case);
s_Lu = zeros(1, n_case);
s_Ld = zeros(1, n_case);
stdat = cell(1, n_case);

for cid = 1:n_case
  case_th = s_case_th(cid);

clear('mode_ST', 'mode_IF', 'pE', 'pI');  % if exist('mode_ST','var')
res_dir = 'result_GC_simple/';
switch case_th
  case 1
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 2
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);
  case 3
    load([res_dir,'result_IF_net_100_01_sc=0.005_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 4
    load([res_dir,'result_IF_net_100_01_sc=0.005_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);
  case 5
    load([res_dir,'result_IF_ST_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 6
    load([res_dir,'result_IF_ST_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);
  case 7
    load([res_dir,'result_IF_ST_net_100_01_sc=0.005_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 8
    load([res_dir,'result_IF_ST_net_100_01_sc=0.005_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);
  case 9
    disp('user define case 1: sync');
    load([res_dir,'result_IF_ST_net_100_01_sc=0.005_pr=6.00_ps=0.001_stv=0.50_t=1.00e+06.mat']);
  case 10
    disp('user define case 2: sync');
    load([res_dir,'result_IF_ST_net_100_01_sc=0.005_pr=6.00_ps=0.001_stv=0.50_t=1.00e+07.mat']);
  case 12
    disp('user define case 4: sync');
    load([res_dir,'result_IF_net_100_01_sc=0.005_pr=6.00_ps=0.001_stv=0.50_t=1.00e+07.mat']);
  case 20
    disp('stv = 0.0625');
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=0.06_t=1.25e+06.mat']);
  case 21
    disp('stv = 1.0');
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=1.00_t=1.00e+07.mat']);
  case 22
    disp('stv = 1.5');
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=1.50_t=1.00e+07.mat']);
  case 23
    disp('stv = 2.5');
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=2.50_t=1.00e+07.mat']);
  case 24
    disp('stv = 3.5');
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=3.50_t=1.00e+07.mat']);
  case 25
    disp('stv = 4.5');
    load([res_dir,'result_IF_net_100_01_sc=0.010_pr=0.24_ps=0.020_stv=4.50_t=1.00e+07.mat']);
  case 32  % uniform distributed scee
    load([res_dir,'result_IF_net_100_02_sc=0.010_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);
  case 33  % exponential distributed scee
    load([res_dir,'result_IF_net_100_03_sc=0.005_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);
  case 34  % normal distributed scee
    load([res_dir,'result_IF_net_100_04_sc=0.005_pr=0.24_ps=0.020_stv=0.50_t=1.00e+07.mat']);

  case 41  % bad example
    load([res_dir,'result_IF_net_100_10_sc=0.005_pr=1.00_ps=0.012_stv=0.50_t=1.00e+07.mat']);
  case 43  % 
    load([res_dir,'result_IF_net_100_10_sc=0.005_pr=1.00_ps=0.012_stv=0.50_t=1.00e+06.mat']);
  case 44  % 
    load([res_dir,'result_IF_net_100_10_sc=0.005_pr=1.00_ps=0.012_stv=7.50_t=1.00e+06.mat']);

  case 50  % small world, E, stv=0.5
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=0.50_t=1.00e+06.mat']);
  case 51  % small world, E, stv=1
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=1.00_t=1.00e+06.mat']);
  case 52  % small world, E, stv=2
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=2.00_t=1.00e+06.mat']);
  case 53  % small world, E, stv=3
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=3.00_t=1.00e+06.mat']);
  case 54  % small world, E, stv=4
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=4.00_t=1.00e+06.mat']);
  case 55  % small world, E, stv=5
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=5.00_t=1.00e+06.mat']);
  case 56  % small world, E, stv=6
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=6.00_t=1.00e+06.mat']);
  case 57  % small world, E, stv=7
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=7.00_t=1.00e+06.mat']);
  case 58  % small world, E, stv=7.5
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=7.50_t=1.00e+06.mat']);
  case 59  % small world, E, stv=8
    load([res_dir,'result_IF_net_100_rs01_sc=0.005_pr=1.00_ps=0.012_stv=8.00_t=1.00e+06.mat']);

  case 60  % found a small network ? EI use stv=0.5
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=0.50_t=1.00e+06.mat']);
  case 61  % found a small network ? EI use stv=1
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+06.mat']);
  case 62  % found a small network ? EI use stv=2
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=2.00_t=1.00e+06.mat']);
  case 63  % found a small network ? EI use stv=3
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=3.00_t=1.00e+06.mat']);
  case 64  % found a small network ? EI use stv=4
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=4.00_t=1.00e+06.mat']);
  case 65  % found a small network ? EI use stv=5
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=5.00_t=1.00e+06.mat']);
  case 66  % found a small network ? EI use stv=6
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=6.00_t=1.00e+06.mat']);
  case 67  % found a small network ? EI use stv=7
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=7.00_t=1.00e+06.mat']);
  case 68  % found a small network ? EI use stv=8
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=8.00_t=1.00e+06.mat']);
  case 69  % found a small network ? EI use stv=9
    load([res_dir,'result_IF_net_100_rs01_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=9.00_t=1.00e+06.mat']);

  % extra cases about E and I
  case 101
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.004,0.004,0.008,0.008]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+06.mat']);
  case 102
    load([res_dir,'result_IF_net_100_21_p[80,20]_sc=[0.004,0.004,0.008,0.008]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+06.mat']);
  case 103
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.004,0.004,0.008,0.008]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+07.mat']);
  case 104
    load([res_dir,'result_IF_net_100_21_p[80,20]_sc=[0.004,0.004,0.008,0.008]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+07.mat']);
  case 105
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=0.50_t=1.00e+06.mat']);
  case 106
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+06.mat']);
  case 107
    load([res_dir,'result_IF_net_100_21_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=0.50_t=1.00e+06.mat']);
  case 108
    load([res_dir,'result_IF_net_100_21_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=1.00_t=1.00e+06.mat']);
  case 109
    load([res_dir,'result_IF_net_100_01_p[80,20]_sc=[0.006,0.006,0.006,0.006]_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 110 % C.4a
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.006,0.006,0.006,0.006]_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 111 % C.5a
    load([res_dir,'result_IF_net_100_21_p[80,20]_sc=[0.006,0.006,0.006,0.006]_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 112 % C.4c
    load([res_dir,'result_IF_ST_net_100_20_p[80,20]_sc=[0.006,0.006,0.006,0.006]_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 113 % C.5c
    load([res_dir,'result_IF_ST_net_100_21_p[80,20]_sc=[0.006,0.006,0.006,0.006]_pr=0.24_ps=0.020_stv=0.50_t=1.00e+06.mat']);
  case 114 % C.2c
    load([res_dir,'result_IF_ST_net_100_21_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=1.00_ps=0.012_stv=0.50_t=1.00e+06.mat']);
  case 115 % C.6
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.004,0.004,0.007,0.007]_pr=1.00_ps=0.005_stv=0.50_t=1.00e+06.mat']);
  case 116 % C.6
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=0.50_ps=0.009_stv=0.50_t=1.00e+06.mat']);
  case 117 % C.6
    load([res_dir,'result_IF_net_100_20_p[80,20]_sc=[0.005,0.005,0.007,0.007]_pr=0.42_ps=0.010_stv=0.50_t=1.00e+06.mat']);

  otherwise
    error('no this case!');
end

neu_network = getnetwork(netstr);
p = size(neu_network,1);

% setup output path and format
pic_prefix = '';
if ~exist('mode_IF','var') || strcmpi(mode_IF,'IF')
    pic_prefix = [pic_prefix0, 'IF'];
else
    pic_prefix = [pic_prefix0, 'EIF'];
end
if exist('mode_ST','var') && mode_ST
    pic_prefix = [pic_prefix, '_ST'];
end
if ~exist('pI', 'var') || pI==0
  st_para = sprintf('_%s_sc=%1.3f_pr=%1.2f_ps=%1.3f_stv=%1.2f_t=%1.2e',...
                    netstr, scee, pr, ps, stv, simu_time);
else
  if exist('sc','var')
    st_para = sprintf('_%s_p[%d,%d]_sc=[%1.3f,%1.3f,%1.3f,%1.3f]_pr=%1.2f_ps=%1.3f_stv=%1.2f_t=%1.2e',...
                      netstr, pE, pI, sc(1), sc(2), sc(3), sc(4), pr, ps, stv, simu_time);
  else
    st_para = sprintf('_%s_p[%d,%d]_sc=%1.3f_pr=%1.2f_ps=%1.3f_stv=%1.2f_t=%1.2e',...
                      netstr, pE, pI, scee, pr, ps, stv, simu_time);
  end
end
pic_prefix = [pic_prefix, st_para];
pic_output       = @(st)print('-deps'  ,[pic_prefix, st, '.eps']);
pic_output_color = @(st)print('-depsc2',[pic_prefix, st, '.eps']);
%pic_output       = @(st)print('-dpng'  ,[pic_prefix, st, '.png']);
%pic_output_color = @(st)print('-dpng',[pic_prefix, st, '.png']);
%pic_output       = @(st) [];
%pic_output_color = @(st) [];
% for save data
pic_data_save = @(st, varargin)save('-v7', [pic_prefix, st, '.mat'], 'varargin');
%pic_data_save = @(st, varargin) [];

disp(pic_prefix);

use_od = bic_od;
gc_zero_line = chi2inv(1-p_val, use_od)/len;
clear('gc_zero_cut','scale_gc','s_gc_tick','s_gc_tick_white','s_gc_hist_tick','s_ISI_tick','s_ISI_tick_white');
%{
switch case_th
  case 1
    %gc_zero_cut = 5.5 * gc_zero_line;        % S=0.010, t=1e6
    scale_gc = 10^-3;
    s_gc_tick = [0 1 2 3];
    s_gc_tick_white = [0.5 1.5 2.5];
    s_gc_hist_tick = [0 40 80 120 160];
    s_ISI_tick = [4 5 6 7];
    s_ISI_tick_white = 0.5+[4 5 6 7];
  case 2
    %gc_zero_cut = 39 * gc_zero_line;         % S=0.010, t=1e7
    scale_gc = 10^-3;
    s_gc_tick = [0 1 2 3];
    s_gc_tick_white = [0.5 1.5 2.5];
    s_gc_hist_tick = [0 40 80 120 160];
    s_ISI_tick = [4 5 6 7];
    s_ISI_tick_white = 0.5+[4 5 6 7];
  case 3
    %gc_zero_cut = 1.015 * gc_zero_line;      % S=0.005, t=1e6
    scale_gc = 10^-5;
    s_gc_tick = [0 2 4 6];
    s_gc_tick_white = [1 3 5];
    s_gc_hist_tick = [0 75 150 225];
    s_ISI_tick = [40 50 60 70];
    s_ISI_tick_white = 5+s_ISI_tick;
  case 4
    %gc_zero_cut = 4 * gc_zero_line;          % S=0.005, t=1e7
    scale_gc = 10^-5;
    s_gc_tick = [0 1 2 3 4];
    s_gc_tick_white = [];
    s_gc_hist_tick = [0 40 80 120];
    s_ISI_tick = [40 50 60 70];
    s_ISI_tick_white = 5+s_ISI_tick;
  case 5
    %gc_zero_cut = 5.0 * gc_zero_line;        % S=0.010, t=1e6, ST
    scale_gc = 10^-3;
    s_gc_tick = [0 1 2];
    s_gc_tick_white = [0.5 1.5];
    s_gc_hist_tick = [0 20 40 60 80];
    s_ISI_tick = [4 5 6 7];
    s_ISI_tick_white = 0.5+[4 5 6 7];
  case 6
    %gc_zero_cut = 42.0 * gc_zero_line;       % S=0.010, t=1e7, ST
    scale_gc = 10^-3;
    s_gc_tick = [0 1 2];
    s_gc_tick_white = [0.5 1.5];
    s_gc_hist_tick = [0 20 40 60 80];
    s_ISI_tick = [4 5 6 7];
    s_ISI_tick_white = 0.5+[4 5 6 7];
  case 7
    %gc_zero_cut = 2.2 * gc_zero_line;        % S=0.005, t=1e6, ST
    scale_gc = 10^-5;
    s_gc_tick = [0 2 4 6];
    s_gc_tick_white = [1 3 5];
    s_gc_hist_tick = [0 20 40 60 80];
    s_ISI_tick = [40 50 60 70];
    s_ISI_tick_white = 5+s_ISI_tick;
  case 8
    %gc_zero_cut = 2 * gc_zero_line;          % S=0.005, t=1e7, ST
    scale_gc = 10^-4;
    s_gc_tick = [0 0.5 1.0 1.5];
    s_gc_tick_white = [];
    s_gc_hist_tick = [0 40 80 120];
    s_ISI_tick = [40 50 60 70];
    s_ISI_tick_white = 5+s_ISI_tick;
  case 9
    %gc_zero_cut = 64 * gc_zero_line;
    scale_gc = 10^-3;
    s_gc_tick = [0 2.0 4.0 6.0];
    s_gc_tick_white = [s_gc_tick(1:end)+1.0];
    s_gc_hist_tick = [0 80 160 240];
    s_ISI_tick = [48 52 56 60];
    s_ISI_tick_white = diff(s_ISI_tick(1:2))/2+s_ISI_tick;
  case 10
    %gc_zero_cut = 402 * gc_zero_line;
    scale_gc = 10^-3;
    s_gc_tick = [0 2.0 4.0 6.0];
    s_gc_tick_white = [s_gc_tick(1:end)+1.0];
    s_gc_hist_tick = [0 80 160 240];
    s_ISI_tick = [48 52 56 60];
    s_ISI_tick_white = diff(s_ISI_tick(1:2))/2+s_ISI_tick;
  otherwise
    %gc_zero_cut = 1.0 * gc_zero_line;
    scale_gc = 10^-3;
end
%}

if b_use_pairGC
  GC = pairRGrangerT(R(:,1:(use_od+1)*p));
else  
  GC = oGC(:,:,use_od);
end

% auto determine cutline
plain_gc = GC;
plain_gc = plain_gc(eye(p)==0);
plain_network = neu_network;
plain_network = plain_network(eye(p)==0);

[sorted_gc, sorted_gc_id] = sort(plain_gc);
sorted_gc_zeroone = zeros(1, length(plain_gc));
sorted_gc_zeroone = plain_network(sorted_gc_id);
tmp_0 = cumsum(sorted_gc_zeroone==0);
tmp_0 = tmp_0(end) - tmp_0;
tmp_1 = cumsum(sorted_gc_zeroone>0);
[min_err_num, min_err_id] = min([tmp_0; 0] + [0; tmp_1]);
min_err_gc = sorted_gc(min_err_id);
%disp(['min_err_gc = ', num2str(min_err_gc)]);
p_val_min_err = 1-chi2cdf(min_err_gc*len, use_od);

if auto_gc_zero_cut
  gc_zero_cut = min_err_gc;
  st_auto_gc1 = ''
  st_auto_gc2 = '  (in use)';
else
  gc_zero_cut = gc_zero_line;
  st_auto_gc1 = '  (in use)';
  st_auto_gc2 = '';
end

if exist('s_gc_tick','var')
  s_gc_tick = s_gc_tick';
end
if exist('s_gc_hist_tick','var')
  s_gc_hist_tick = s_gc_hist_tick';
end
if exist('s_gc_tick_white','var')
  s_gc_tick_white = s_gc_tick_white';
end
if exist('s_ISI_tick','var')
  s_ISI_tick = s_ISI_tick';
end
if exist('s_ISI_tick_white','var')
  s_ISI_tick_white = s_ISI_tick_white';
end
%neu_network_guess = gc_prob_nonzero(GC,use_od,len) >= 1 - p_val;
neu_network_guess = GC >= gc_zero_cut;
adj_cmp = neu_network_guess - neu_network;

disp(['aic_od = ',num2str(aic_od)]);
disp(['bic_od = ',num2str(bic_od)]);
disp(['use_od = ',num2str(use_od)]);
mean_gc0 = mean(GC(neu_network==0));
var_gc0  = var (GC(neu_network==0));
var_theo_gc0 = 2*use_od/len^2 + 4/len*(mean_gc0-use_od/len);
mean_gc1 = mean(GC(neu_network~=0));
var_gc1  = var (GC(neu_network~=0));
var_theo_gc1 = 2*use_od/len^2 + 4/len*(mean_gc1-use_od/len);

fprintf('data length: %.2e\n', len);
scale_gc = 1e-5;
if scale_gc~=0
fprintf('mean GC0/%.1e: %.2f  (stdv=%.2f theo: %.2f res stdv=%.2f)\n',...
        scale_gc, mean_gc0/scale_gc, sqrt(var_gc0)/scale_gc,...
        sqrt(var_theo_gc0)/scale_gc, sqrt(var_gc0-var_theo_gc0)/scale_gc);
fprintf('mean GC1/%.1e: %.2f  (stdv=%.2f theo: %.2f res stdv=%.2f)\n',...
        scale_gc, mean_gc1/scale_gc, sqrt(var_gc1)/scale_gc,...
        sqrt(var_theo_gc1)/scale_gc, sqrt(var_gc1-var_theo_gc1)/scale_gc);

stdat{cid} = sprintf('%.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n',...
        use_od/len/scale_gc, 
        mean_gc0/scale_gc, mean_gc1/scale_gc, sqrt(var_gc0)/scale_gc,...
        sqrt(var_theo_gc0)/scale_gc, sqrt(var_gc1)/scale_gc,...
        sqrt(var_theo_gc1)/scale_gc);
else
fprintf('mean GC0: %.2e  (stdvar=%.2e, theo: %.2e)  res stdvar=%.2e\n',...
        mean_gc0, sqrt(var_gc0), sqrt(var_theo_gc0),...
        sqrt(var_gc0-var_theo_gc0));
fprintf('mean GC1: %.2e  (stdvar=%.2e, theo: %.2e)  res stdvar=%.2e\n',...
        mean_gc1, sqrt(var_gc1), sqrt(var_theo_gc1),...
        sqrt(var_gc1-var_theo_gc1));
end
fprintf('   ratio: %.2f\n', mean_gc1/mean_gc0);
fprintf('Ideal GC0 = %.2e\n', use_od/len);
fprintf('GC(P=%.2e) = %.2e%s\n',...
        p_val,         gc_zero_line, st_auto_gc1);
fprintf('GC(P=%.2e) = %.2e  (least wrong edges)%s\n',...
        p_val_min_err, min_err_gc,   st_auto_gc2);

n_over = sum((adj_cmp(:)==1));
n_lack = sum((adj_cmp(:)==-1));
n_conn = sum(neu_network(:));
disp(['over guess: ', num2str(n_over)]);
disp(['lack guess: ', num2str(n_lack)]);
disp(['total conn: ', num2str(n_conn), '  max: ',num2str(p*(p-1))]);

if b_cal_net
    addpath('/home/xyy/matcode/octave-networks-toolbox');
    addpath('/home/xyy/matcode/prj_neu_smallworld');
    s_Ld(cid) = avePathLength(neu_network_guess);
    s_Lu(cid) = avePathLength((neu_network_guess+neu_network_guess.')>0);
    s_CC(cid) = clustering_coef_BDN(neu_network_guess);
end
s_over(cid) = n_over;
s_lack(cid) = n_lack;
s_stv(cid)  = stv;
s_aic(cid)  = aic_od;
s_bic(cid)  = bic_od;

if b_output_pics

% show histogram of ISI
figure(1);  set(gca,'fontsize',font_size);
hist(ISI, 20);
ylabel('count');
xlabel('ave ISI of each neuron (ms)');
if exist('s_ISI_tick','var') && ~isempty(s_ISI_tick)
    set(gca,'xtick',[s_ISI_tick; s_ISI_tick_white],'xticklabel',[num2str(s_ISI_tick); char(32*ones(length(s_ISI_tick_white),1))]);
end
pic_output_color('_ISI_hist');

% show true adjacency matrix
[xx,yy] = meshgrid(0:p, 0:p);
figure(2);  set(gca,'fontsize',font_size);
colormap(fliplr(gray()')');
shading('faceted');
%hd=pcolor(xx,yy,[[flipud(neu_network); zeros(1,p)],zeros(p+1,1)]);
hd=pcolor(xx,yy,[[neu_network; zeros(1,p)],zeros(p+1,1)]);
hd=colorbar();
set(hd,'fontsize',font_size);
set(hd,'ytick',[0;0.5;1; 0.25;0.75],'yticklabel',['0';'0.5';'1.0';' ';' ']);
%axis('square','off');
xlabel('index of driven one');
ylabel('index of passive one');
set(gca,'xtick',[1 50 100]);
set(gca,'ytick',[1 50 100]);
pic_output('_gcadj_ans');

% show GC adjacency matrix (gray points)
figure(3);  set(gca,'fontsize',font_size);
colormap(fliplr(gray()')');
shading('faceted');
adj_gc = zero_GC/scale_gc;
adj_gc(adj_gc<0) = 0;
%hd=pcolor(xx,yy,[[flipud(adj_gc); zeros(1,p)],zeros(p+1,1)]);
hd=pcolor(xx,yy,[[adj_gc; zeros(1,p)],zeros(p+1,1)]);
%caxis([0,2]);
hd=colorbar();
%set(hd,'ylabel','GC');
set(hd,'fontsize',font_size);
%set(hd,'ytick',s_gc_tick);
if exist('s_gc_tick','var') && ~isempty(s_gc_tick)
  set(hd,'ytick',[s_gc_tick; s_gc_tick_white],'yticklabel',[num2str(s_gc_tick); char(32*ones(length(s_gc_tick_white),1))]);
end
%axis('square','off');
xlabel('index of driven one');
ylabel('index of passive one');
set(gca,'xtick',[1 50 100]);
set(gca,'ytick',[1 50 100]);
%zlabel('GC value (10^{',num2str(round(log10(scale_gc))),'})');
pic_output('_gcadj');

% show correctness by color map
figure(4);  set(gca,'fontsize',font_size);
shading('faceted');
%pcolor(xx,yy,[[flipud(adj_cmp); zeros(1,p)],zeros(p+1,1)]);
pcolor(xx,yy,[[adj_cmp; zeros(1,p)],zeros(p+1,1)]);
caxis([-1,1]);
hd=colorbar();
set(hd,'fontsize',font_size);
%axis('square','off');
xlabel('index of driven one');
ylabel('index of passive one');
set(gca,'xtick',[1 50 100]);
set(gca,'ytick',[1 50 100]);
pic_output_color('_adj_cmp');

% show sorted GC values
figure(5);  set(gca,'fontsize',font_size);
hd=plot([1,(p*p-p)]/10^3, [gc_zero_line, gc_zero_line]/scale_gc, 'b-',...
        [1,(p*p-p)]/10^3, [  min_err_gc,   min_err_gc]/scale_gc, 'r--',...
        [1:(p*p-p)]/10^3, sort(plain_gc)/scale_gc, 'k.');
set(hd,'linewidth',line_width);
xlabel('sorted index (10^3)');
ylabel(['GC value (10^{',num2str(round(log10(scale_gc))),'})']);
%set(gca,'ytick',s_gc_tick);
if exist('s_gc_tick') && ~isempty('s_gc_tick')
    set(gca,'ytick',[s_gc_tick; s_gc_tick_white],'yticklabel',[num2str(s_gc_tick); char(32*ones(length(s_gc_tick_white),1))]);
end
hd=legend(['GC0\neq0 P-value = ',num2str(p_val)], 'best 0 or 1 threshold',...
          'location','northwest');
set(hd,'fontsize',font_size-6);
pic_output_color('_gc_sort');
pic_data_save('_gc_sort', p, use_od, len, plain_gc, plain_network);

% show histogram of GC values
figure(6);  set(gca,'fontsize',font_size);
[nn,xx] = hist(plain_gc/scale_gc,100);
peak_nn = max(nn(10:end));
line_high = peak_nn;
peak_scale = 10^floor(log10(1.3*peak_nn));
peak_scale = ceil(4*1.3*peak_nn/peak_scale)/4*peak_scale;
nn_id = 1;
while nn(nn_id)>peak_scale
    disp(['bin[',num2str(nn_id),']=',num2str(nn(nn_id))]);
    line_high = peak_scale;
    nn_id = nn_id + 1;
end

bar(xx,nn, 'facecolor',[0.5,0.5,0.5],'linestyle','-');
xlabel(['GC value (10^{',num2str(round(log10(scale_gc))),'})']);
ylabel('count/bin');
%set(gca, 'xtick',s_gc_tick);
if exist('s_gc_tick') && ~isempty('s_gc_tick')
    set(gca,'xtick',[s_gc_tick; s_gc_tick_white],'xticklabel',[num2str(s_gc_tick); char(32*ones(length(s_gc_tick_white),1))]);
end
if exist('s_gc_hist_tick') && ~isempty('s_gc_hist_tick')
    set(gca,'ytick',s_gc_hist_tick);
end
hold('on');
hd=plot([gc_zero_line,gc_zero_line]/scale_gc, [0,line_high],'g-',...
        [  min_err_gc,  min_err_gc]/scale_gc, [0,line_high],'r--');
if line_high == peak_scale
    sa=axis();  sa(4)=line_high;  axis(sa);
end
set(hd,'linewidth',line_width);
hd=legend('', ['GC0\neq0 P-value = ',num2str(p_val)], 'best 0 or 1 threshold', 'asdf');
set(hd,'fontsize',font_size-6);
hold('off');
pic_output_color('_gc_hist');

end  % if b_output_pics
%test_big_p

end %for case_th

for id=1:length(stdat)
  fprintf('%.1f  %d(%d)  %s', s_stv(id), s_bic(id), s_aic(id), stdat{id});
end

for id=1:length(s_CC)
  fprintf('%.1f  %d(%d)  %d  %d  %.2f  %.2f  %.3f\n',...
          s_stv(id), s_bic(id), s_aic(id), s_over(id), s_lack(id), s_Ld(id), s_Lu(id), s_CC(id));
end

