% show scan: percentage of common poisson input, fix others
tic;

%indir  = 'data_scan/';
%fname  = 'comm_w01_test_net_2_2_sc=0.01_pr=1.00_ps=0.01_t=1.0e+05_stv=5.0e-01.mat';

indir  = 'data_scan_comm/';
%fname  = 'comm_w01_net_2_2_sc=0.01_pr=1.00_ps=0.01_t=1.0e+06_stv=5.0e-01.mat';
%fname  = 'comm_w02_net_2_2_sc=0.01_pr=1.00_ps=0.01_t=2.0e+06_stv=1.0e+00.mat';
%fname  = 'comm_w02_net_2_2_sc=0.01_pr=1.00_ps=0.01_t=1.0e+07_stv=5.0e-01.mat';

%net_2_2, 0.01, 1.0, 0.010, 0.5, 1.0/32, 1, 1e7  w02  12567.7sec @xhp

s_file = {
'comm_w11_net_2_2_sc=0.01_pr=8.00_ps=0.005_t=1.0e+07_stv=5.0e-01.mat',
'comm_w12_net_2_2_sc=0.01_pr=4.00_ps=0.01_t=1.0e+07_stv=5.0e-01.mat',
'comm_w13_net_2_2_sc=0.01_pr=2.00_ps=0.02_t=1.0e+07_stv=5.0e-01.mat',
'comm_w14_net_2_2_sc=0.01_pr=2.00_ps=0.005_t=1.0e+07_stv=5.0e-01.mat',
'comm_w15_net_2_2_sc=0.01_pr=1.00_ps=0.01_t=1.0e+07_stv=5.0e-01.mat',
'comm_w16_net_2_2_sc=0.01_pr=0.50_ps=0.02_t=1.0e+07_stv=5.0e-01.mat',
'comm_w17_net_2_2_sc=0.01_pr=1.00_ps=0.005_t=1.0e+07_stv=5.0e-01.mat',
'comm_w18_net_2_2_sc=0.01_pr=0.50_ps=0.01_t=1.0e+07_stv=5.0e-01.mat',
'comm_w19_net_2_2_sc=0.01_pr=0.25_ps=0.02_t=1.0e+07_stv=5.0e-01.mat',

'comm_w21_ST_net_2_2_sc=0.01_pr=8.00_ps=0.005_t=1.0e+07_stv=5.0e-01.mat',
'comm_w22_ST_net_2_2_sc=0.01_pr=4.00_ps=0.01_t=1.0e+07_stv=5.0e-01.mat',
'comm_w23_ST_net_2_2_sc=0.01_pr=2.00_ps=0.02_t=1.0e+07_stv=5.0e-01.mat',
'comm_w24_ST_net_2_2_sc=0.01_pr=2.00_ps=0.005_t=1.0e+07_stv=5.0e-01.mat',
'comm_w25_ST_net_2_2_sc=0.01_pr=1.00_ps=0.01_t=1.0e+07_stv=5.0e-01.mat',
'comm_w26_ST_net_2_2_sc=0.01_pr=0.50_ps=0.02_t=1.0e+07_stv=5.0e-01.mat',
'comm_w27_ST_net_2_2_sc=0.01_pr=1.00_ps=0.005_t=1.0e+07_stv=5.0e-01.mat',
'comm_w28_ST_net_2_2_sc=0.01_pr=0.50_ps=0.01_t=1.0e+07_stv=5.0e-01.mat',
'comm_w29_ST_net_2_2_sc=0.01_pr=0.25_ps=0.02_t=1.0e+07_stv=5.0e-01.mat'
};

%net_2_2, 0.01, 8.0, 0.005, 0.5, 1.0/32, 1, 1e7  w11  13658.6sec @xhp
%net_2_2, 0.01, 4.0, 0.010, 0.5, 1.0/32, 1, 1e7  w12  12745.6sec @xhp
%net_2_2, 0.01, 2.0, 0.020, 0.5, 1.0/32, 1, 1e7  w13  12235.3sec @xhp
%net_2_2, 0.01, 2.0, 0.005, 0.5, 1.0/32, 1, 1e7  w14  12666.6sec @dell
%net_2_2, 0.01, 1.0, 0.010, 0.5, 1.0/32, 1, 1e7  w15  12340.0sec @dell
%net_2_2, 0.01, 0.5, 0.020, 0.5, 1.0/32, 1, 1e7  w16  12235.3sec @dell
%net_2_2, 0.01, 1.0, 0.005, 0.5, 1.0/32, 1, 1e7  w17  12580.6sec @dell  "The neuron does not fire at the end point?!\nvoltage difference at the beginning time and ending time: 1.33893e-13;4.53264e-07"
%net_2_2, 0.01, 0.5, 0.010, 0.5, 1.0/32, 1, 1e7  w18  12454.7sec @dell
%net_2_2, 0.01, 0.25,0.020, 0.5, 1.0/32, 1, 1e7  w19  12367.6sec @dell


pic_prefix0 = 'pic_tmp/';

if ~exist('font_size','var')
  font_size = 24;
end
if ~exist('line_width','var')
  line_width = 3;
end

od_mode = 1;               % 1 is 'BIC', 2 is 'AIC', 3 is 'BICall'
switch od_mode
  case 1
    f_od_mode = @(vBIC,vAIC,vBICall) vBIC;
  case 2
    f_od_mode = @(vBIC,vAIC,vBICall) vAIC;
  case 3
    f_od_mode = @(vBIC,vAIC,vBICall) vBICall;
  otherwise
    error('invalid od mode');
end

for id_file = 1:length(s_file);
    fname = s_file{id_file};
    load([indir,fname]);

    % ignore some of them
%    s_common_cut = 0.7;
    s_common_cut = 1.0;
    s_ignore = s_common_poisson>s_common_cut;
    s_common_poisson(s_ignore) = [];
    s_aveISI (:,  s_ignore) = [];
    s_ISI_dis(:,:,s_ignore) = [];

    %m_aic = zeros(1,length(s_common_poisson));
    %m_bic = zeros(1,length(s_common_poisson));
    %m_bic_all = zeros(1,length(s_common_poisson));
    m_GC  = zeros(p,p,length(s_common_poisson));
    m_lGC = zeros(p,p,length(s_common_poisson));
    m_uGC = zeros(p,p,length(s_common_poisson));
    m_De  = zeros(p,p,length(s_common_poisson));
    m_inst_cor = zeros(1, length(s_common_poisson));
    for id_s_comm = 1:length(s_common_poisson)
        oGC     = s_oGC    (:,:,:,id_s_comm);
        oDe     = s_oDe    (:,:,:,id_s_comm);
        [aic_od, bic_od, zero_GC, oAIC, oBIC] = AnalyseSeries2(s_od, oGC, oDe, len);
    %    R       = s_R      (:,:,  id_s_comm);
    %    [od_joint, od_vec] = chooseROrderFull(R, len, 'BIC');
    %    bic_od_all = max([od_joint, od_vec]);
    %    m_aic(id_s_comm)     = aic_od;
    %    m_bic(id_s_comm)     = bic_od;
    %    m_bic_all(id_s_comm) = bic_od_all;
        bic_od_all = 0;  % doesn't use it

        od = f_od_mode(bic_od, aic_od, bic_od_all);  % bic_od;
        GC = oGC(:,:,od);                            % or zero_GC;  oGC(:,:,20);  oGC(:,:,bic_od);
        [lGC, uGC] = gc_prob_intv(GC, od, len);

        m_GC (:,:,id_s_comm) = GC;
        m_lGC(:,:,id_s_comm) = lGC;  % lower bound
        m_uGC(:,:,id_s_comm) = uGC;  % upper bound

        De = oDe(:,:,od);
        m_De (:,:,id_s_comm) = De;
        m_inst_cor(id_s_comm)= De(1,2)/sqrt(De(1,1)*De(2,2));  % instantaneous correlation
    %    aveISI  = s_aveISI (:,    id_s_comm);
    %    ISI_dis = s_ISI_dis(:,:,  id_s_comm);
    end

fprintf('net:%s, sc:%.3f, pr:%.2f, ps:%.4f, time:%.2e,stv:%.2f,len:%.2e\n',...
        netstr, scee, pr, ps, simu_time, stv, len);
disp(['ave ISI: ',num2str(mean(s_aveISI(1,:)))]);
%    toc;

    % graph output functions
    if mode_IF
        pic_prefix = [pic_prefix0, 'IF'];
    else
        pic_prefix = [pic_prefix0, 'expIF'];
    end
    if mode_ST
        pic_prefix = [pic_prefix, '_ST'];
    end
    pic_prefix = sprintf('%s_%s_sc=%.3f_pr=%.2f_ps=%.3f_', pic_prefix, netstr, scee, pr, ps);
    pic_suffix = sprintf('_fullcommon_stv=%.2f_t=%.2e', stv, simu_time);
    %pic_output = @(st)print('-dpng',[pic_prefix, st, pic_suffix, '.png'],'-r100');    % output function
    %pic_output_color = pic_output;                                        % for color output
    pic_output       = @(st)print('-deps'  ,[pic_prefix, st, pic_suffix, '.eps']);
    pic_output_color = @(st)print('-depsc2',[pic_prefix, st, pic_suffix, '.eps']);

    % Generate labels in legend()
    label_st_ISI = {};
    label_st_freq = {};
    for kk=1:p
        label_st_ISI  = {label_st_ISI{:}, ['ISI ',num2str(kk)]};
        label_st_freq = {label_st_freq{:}, ['Freq ',num2str(kk)]};
    end

    % ISI v.s.  % of common poisson input
    figure(1);  set(gca, 'fontsize',font_size);
    plot(100.0*s_common_poisson, s_aveISI, 'linewidth',line_width);
    hd=legend(label_st_ISI{:});  set(hd, 'fontsize',font_size-2);
    axis([0, 100.0*s_common_cut]);
    xlabel('% of common poisson input');
    ylabel('average ISI/ms');
    pic_output_color('ISI');

    figure(2);  set(gca, 'fontsize',font_size);
    plot(100.0*s_common_poisson, m_inst_cor, '-+', 'linewidth',line_width);
    axis([0, 100.0*s_common_cut]);
    xlabel('% of common poisson input');
    ylabel('cor(\epsilon_t, \eta_t)');
    pic_output_color('inst_cor');

    % convert to 2-D matrix
    plain_GC = reshape(m_GC,p*p,size(m_GC,3));
    plain_GC(1:p+1:end,:) = [];
    plain_lGC = reshape(m_lGC,p*p,size(m_GC,3));
    plain_lGC(1:p+1:end,:) = [];
    plain_uGC = reshape(m_uGC,p*p,size(m_GC,3));
    plain_uGC(1:p+1:end,:) = [];
    plain_lGC = plain_GC - plain_lGC;
    plain_uGC = plain_uGC - plain_GC;

    scale_GC = 1e-3;

    % GC(x->y) & GC(y->x)  v.s.  % of common poisson input
    figure(3);  set(gca, 'fontsize',font_size);
    % get a color map for different curve in plot/errorbar
    %plot(100.0*s_common_poisson, plain_GC/scale_GC);
    hd = errorbar(100.0*s_common_poisson'*ones(1,p*p-p), plain_GC'/scale_GC, plain_lGC'/scale_GC, plain_uGC'/scale_GC);
    if (p>7)
        cm = hsv(p+1);
    else
        cm = [0 0 1; 0 0.5 0; 1 0 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0; 0.25 0.25 0.25];
    end
    for k=1:p
        set(hd(k),'color',cm(k,:), 'linewidth',line_width);
    end
    hd=legend('x->y','y->x');  set(hd, 'fontsize',font_size-2);
    axis([0, 100.0*s_common_cut]);
    xlabel('% of common poisson input');
    ylabel(['GC/',num2str(scale_GC)]);
    pic_output_color('GC_%common');

    % GC(x->y) / GC(y->x)  v.s.  % of common poisson input
    figure(4);  set(gca, 'fontsize',font_size);
    plot(100.0*s_common_poisson, plain_GC(1,:)./plain_GC(2,:), '-+', 'linewidth',line_width);
    hd=legend('GC(x->y) / GC(y->x)');  set(hd, 'fontsize',font_size-2);
    axis([0, 100.0*s_common_cut]);
    xlabel('% of common poisson input');
    ylabel(['GC ratio']);
    pic_output_color('GCratio_%common');

%    toc;
end

