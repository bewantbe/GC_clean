% 
%addpath('../');
%addpath('../../GCcal');
use_od = 17;

mode_IF = 'IF';
mode_ST = 0;
netstr = 'net_2_2';
scee = 0.01;
pr = 1;
ps = 0.012;
simu_time = 1e5;
extst = 'new --RC-filter -q';

if ~exist('extst','var')
    extst = '';
end
if strcmpi(mode_IF,'ExpIF')
    extst = ['ExpIF ',extst];
end
if exist('dt','var')
    extst = [extst, sprintf(' -dt %.16e',dt)];
end
if ~exist('stv','var')
    stv = 1/2;
end

tic
n_stat = 14000;
s_gc_stat = zeros(n_stat, 2);
for k=1:n_stat
    init_len = 1e4;
    X = gendata_neu(netstr, scee, pr, ps, simu_time+init_len*stv, stv, extst);
    X = X(:, init_len+1:end);
    gc = nGrangerT(X, use_od);
    s_gc_stat(k,:) = [gc(2), gc(3)];
    fprintf('k=%d\n', k);  fflush(stdout);
end
save('-v7', 's_gc_stat.mat', 's_gc_stat');
len = round(simu_time/stv);
save('-v7', sprintf('s_gc_stat=%d_od=%d_L=%.1e_pr=%.1e_ps=%.1e_sc=%.1e.mat',...
     n_stat, use_od, len, pr, ps, scee),...
     's_gc_stat','len','use_od','pr','ps','scee','netstr','stv',...
     'mode_IF','mode_ST','simu_time','extst');
toc


%s_gc_stat_w4 = s_gc_stat;
%save('-v7', 's_gc_stat_w4.mat', 's_gc_stat_w4');

%load w1/s_gc_stat_w1.mat
%load w2/s_gc_stat_w2.mat
%load w3/s_gc_stat_w3.mat
%load w4/s_gc_stat_w4.mat
%s_gc_stat = cat(1, s_gc_stat_w1, s_gc_stat_w2, s_gc_stat_w3, s_gc_stat_w4);

%len = round(simu_time/stv);
%save(sprintf('s_gc_stat_od=%d_L=%.1e_pr=%.1e_ps=%.1e_sc=%.1e.mat',...
%     use_od, len, pr, ps, scee),...
%     's_gc_stat','len','use_od','pr','ps','scee','netstr','stv',...
%     'mode_IF','mode_ST','simu_time','extst');

%figure(1);
%gc_scale = 1e-4;
%n_bin = 50;
%[nn1, xx1] = hist_pdf(s_gc_stat(:,1)/gc_scale, n_bin);
%[nn0, xx0] = hist_pdf(s_gc_stat(:,2)/gc_scale, n_bin);
%stairs([0, xx1], [0, nn1], 'color', 'blue');
%hold on
%stairs([0, xx0], [0, nn0], 'color', [0 0.7 0]);
%hold off

%load('s_gc_stat_od=17_L=2.0e+05_pr=1.0e+00_ps=1.2e-02_sc=1.0e-02.mat');
%gc_scale = 1e-4;
%g_max    = 5;
%nbins    = 40;
%figure(1);
%nn = hist2dnn(s_gc_stat(:,2)/gc_scale, s_gc_stat(:,1)/gc_scale, [0 g_max], nbins);
%rg0 = g_max*(1:nbins+1)/nbins;
%rg1 = g_max*(1:nbins+1)/nbins;
%pcolor(rg0, rg1, [[nn, zeros(length(rg1)-1,1)]; zeros(1,length(rg0))]);
%colormap(flipud(colormap('gray')));
%shading('interp');
%axis([0 g_max 0 g_max]);

%figure(2);
%s_g = g_max*(1:nbins)/nbins;
%nn0 = sum(nn, 2)';
%nn1 = sum(nn, 1);
%nn01 = nn0'*nn1;
%rg0 = g_max*(1:nbins+1)/nbins;
%rg1 = g_max*(1:nbins+1)/nbins;
%pcolor(rg0, rg1, [[nn01, zeros(length(rg1)-1,1)]; zeros(1,length(rg0))]);
%colormap(flipud(colormap('gray')));
%shading('interp');
%axis([0 g_max 0 g_max]);

