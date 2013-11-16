% copy from Chap 3. Several Examples
% pic ratio of minimum causality
tic();

% for final thesis
s_signature = {'data_scan_all_net_2_2_t1e6/sc_3'};
%s_signature = {'data_scan_all_net_3/w_04_07_sc2'};

%s_signature = {'data_scan_all_net_3/w_00_03', 'data_scan_all_net_3/w_04_07', 'data_scan_all_net_3/w_08_11', 'data_scan_all_net_3/w_12_15'};
%s_signature = {'data_scan_all_net_3/w_04_07'};
%s_signature = {'data_scan_all_net_3/w_04_07_sc2'};
%s_signature = {'data_scan_expIF/w1'};
%s_signature = {'data_scan_expIF/w_net_3_06_sc1_t1e6'};  % scan_all_expIF_w106.m
%s_signature = {'data_scan_expIF/w_net_3_06_sc2_t1e6'};  % scan_all_expIF_w126.m
%s_signature = {'data_scan_expIF/w_net_3_06_sc1_t1e5'};  % scan_all_expIF_w86.m
%s_signature = {'data_scan_expIF/w_net_3_06_sc2_t1e5'};  % scan_all_expIF_w66.m
%s_signature = {'data_scan_expIF/w_net_2_2_sc2_t1e6'};   % scan_all_expIF_w200.m  cal in xhp, cost 185715sec

%s_signature = {'data_scan_expIF/w_net_2_2_sc2_t1e6_499'};  %scan_all_expIF_w300.m  cal in dell, cost 259945.9sec
%s_signature = {'data_scan_expIF/w_net_3_09_sc2_t1e6_499'};  %scan_all_expIF_w490.m  cal in dell, cost 306950sec
%s_signature = {'data_scan_expIF/w_net_3_06_sc2_t1e6_499'};  %scan_all_expIF_w226.m  cal in dell, cost 382889.7sec

%s_signature = {'data_scan_expIF/w_net_3_06_sc2_t1e6_499_st'};  %scan_all_expIF_w1460.m  cal in xhp, cost 303532sec  use spike train

%s_signature = {'data_scan/template_ST_expIF_w01'};

%s_signature = {'data_scan_IF/w_net_2_2_sc2_t1e6_199_st'};   %                 3xxxx.x
%s_signature = {'data_scan_IF/w_net_3_06_sc2_t1e6_199_st'};   % Elapsed time is 54839.4 seconds.

pic_prefix0 = 'pic_tmp/';
gc_od_mode = 'BIC';         % the order used for final GC value, 'BIC','AIC','maxBIC','zero'

if isnumeric(gc_od_mode)
  if gc_od_mode==-1
    f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) abs(zero_GC);
  elseif gc_od_mode==0
    f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) zero_GC;
  elseif gc_od_mode>0
    f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) oGC(:,:,gc_od_mode);
  else
    error('invalid GC od mode');
  end
else
  switch gc_od_mode
    case 'BIC'
      f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) oGC(:,:,bic_od);
    case 'AIC'
      f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) oGC(:,:,aic_od);
    case 'maxBIC'
      f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) oGC(:,:,bic_odmax);
    case 'zero'
      f_gc_od = @(oGC, zero_GC, bic_od, aic_od, bic_odmax) zero_GC;
    otherwise
      error('invalid GC od mode');
  end
end
if ~exist('font_size','var')
  font_size = 28;
end
if ~exist('line_width','var')
  line_width = 3;
end

for id_signature = 1:length(s_signature)
signature = s_signature{id_signature};     % to distinguish different parallel program instances (also dir)
% load these variables: 's_net', 's_time', 's_scee', 's_prps', 's_ps', 's_stv', 's_od', 'hist_div', 'maxod'
load([signature, '_info.mat']);
% get the true "save time interval"
if isempty(strfind(signature,'expIF'))
  time_step = 1.0/32;
  p_value = 1e-15;
else
  time_step = 0.004;
  p_value = 1e-15;
end
s_stv = ceil(s_stv/time_step)*time_step;

% show range
s_id_net  = 1:length(s_net);
s_id_time = 1:length(s_time);
s_id_scee = 1:length(s_scee);
s_id_stv  = 1:length(s_stv);
s_id_prps = 1:length(s_prps);
s_id_ps   = 1:length(s_ps);

s_id_stv  = 3;  % overwrite above value
s_id_net  = 1;

for id_net = s_id_net
 netstr = s_net{id_net};
 matname = ['network/', netstr, '.txt'];
 neu_network = load('-ascii', matname);
 p = size(neu_network, 1);
for id_time = s_id_time
 simu_time = s_time(id_time);
for id_scee = s_id_scee
 scee = s_scee(id_scee);

% load these variables: 'prps_ps_stv_oGC', 'prps_ps_stvsignature_oDe', 'prps_ps_stv_R', 'prps_ps_aveISI', 'prps_ps_ISI_dis'
 datamatname = sprintf('%s_%s_sc=%g_t=%.3e.mat', signature, netstr, scee, simu_time);
 load(datamatname);
for id_stv = s_id_stv
 stv = s_stv(id_stv);
 len = round(simu_time/stv);                  % !! I don't know the exact expression

    s_aic_od = zeros(length(s_id_ps), length(s_id_prps));
    s_bic_od = zeros(length(s_id_ps), length(s_id_prps));
    s_all_od = zeros(length(s_id_ps), length(s_id_prps));
    s_zero_GC = zeros(p,p,length(s_id_ps), length(s_id_prps));
    ISI_a_b  = zeros(length(s_id_ps), length(s_id_prps));
    wrong_num = zeros(p,p);

    id_id_prps = 0;
    for id_prps = s_id_prps
     prps = s_prps(id_prps);
     id_id_prps = id_id_prps + 1;
     id_id_ps = 0;
    for id_ps = s_id_ps
     ps = s_ps(id_ps);
     id_id_ps = id_id_ps + 1;
     pr = prps / ps;
        aveISI = prps_ps_aveISI(:, id_prps, id_ps);
        ISI_dis = prps_ps_ISI_dis(:,:,id_prps, id_ps);
        oGC = prps_ps_stv_oGC(:,:,:, id_prps, id_ps, id_stv);
        oDe = prps_ps_stv_oDe(:,:,:, id_prps, id_ps, id_stv);
        [aic_od, bic_od, zero_GC, oAIC, oBIC] = AnalyseSeries2(s_od, oGC, oDe, len);
%        R = prps_ps_stv_R(:,:, id_prps, id_ps, id_stv);
%        [od_joint, od_vec] = chooseROrderFull(R, len, 'BIC');
%        bic_od_all = max([od_joint, od_vec]);
        bic_od_all = 0;
        s_aic_od(id_id_ps, id_id_prps) = aic_od;
        s_bic_od(id_id_ps, id_id_prps) = bic_od;
        s_zero_GC(:,:,id_id_ps, id_id_prps) = f_gc_od(oGC, zero_GC, bic_od, aic_od, bic_od_all);
        ISI_a_b(id_id_ps, id_id_prps) = mean(aveISI);
        net_diff = (gc_prob_nonzero(oGC(:,:,bic_od), bic_od, len)>1-p_value) - neu_network;
        wrong_num(id_id_ps, id_id_prps) = sum(abs(net_diff(:)));
    end  % ps
    end  % prps

    [xx,yy] = meshgrid(s_prps(s_id_prps)/0.001, s_ps(s_id_ps)/0.001);

    s_zero_GC = permute(s_zero_GC, [3,4,1,2]);

    k1 = zeros(size(s_zero_GC,1), size(s_zero_GC,2));
    k2 = zeros(size(s_zero_GC,1), size(s_zero_GC,2));
    k3 = zeros(size(s_zero_GC,1), size(s_zero_GC,2));
    for j1=1:size(s_zero_GC,1)
    for j2=1:size(s_zero_GC,2)
      G1 = s_zero_GC(j1,j2,:,:);
      a = G1(neu_network > 0.5 & eye(p) < 0.5);       % the GC of "one"
      b = G1(neu_network < 0.5 & eye(p) < 0.5);       % the GC of "zero"
      if (isempty(a)==1)
        a = ones(p*(p-1),1)*NaN;
      end
      if (isempty(b)==1)
        b = ones(p*(p-1),1)*NaN;
      end
      amax = max(a);
      amin = min(a);
      bmax = max(b);
      bmin = min(b);
      k1(j1,j2) = amin/bmax;
      k2(j1,j2) = amax/amin;
      k3(j1,j2) = bmax/bmin;
    end
    end
    if a(1)~=a(1) || b(1)~=b(1)
      k1=zeros(size(s_zero_GC,1),size(s_zero_GC,2));
      k2=zeros(size(s_zero_GC,1),size(s_zero_GC,2));
      k3=zeros(size(s_zero_GC,1),size(s_zero_GC,2));
    end

    k1(imag(k1)~=0) = NaN;

    mode_eif= ~isempty(strfind(lower(signature),lower('expIF')));
    if mode_eif
        pic_prefix = [pic_prefix0, 'expIF'];
    else
        pic_prefix = [pic_prefix0, 'IF'];
    end
    mode_st = ~isempty(strfind(lower(signature),lower('ST')));
    if mode_st
        pic_prefix = [pic_prefix, '_ST'];
    end
    pic_prefix = sprintf('%s_%s_sc=%.4f_', pic_prefix, netstr, scee);
    pic_suffix = sprintf('_stv=%.2f_t=%.2e', stv, simu_time);
    %pic_output = @(st)print('-dpng',[pic_prefix, st, pic_suffix, '.png'],'-r100');    % output function
    %pic_output_color = pic_output;                                        % for color output
%    pic_output       = @(st)print('-deps'  ,[pic_prefix, st, pic_suffix, '.eps'], ['-Fontsize:',num2str(font_size)]);
%    pic_output_color = @(st)print('-depsc2',[pic_prefix, st, pic_suffix, '.eps'], ['-Fontsize:',num2str(font_size)]);
    pic_output       = @(st)print('-deps'  ,[pic_prefix, st, pic_suffix, '.eps']);
    pic_output_color = @(st)print('-depsc2',[pic_prefix, st, pic_suffix, '.eps']);

    pic_data_save = @(st, varargin)save('-v7', [pic_prefix, st, '.mat'], 'varargin');

%    figure(1);
%    set(gca, 'fontsize',font_size);
%    pcolor(xx,yy,k1);
%    caxis([0,5]);
%    shading('flat');
%    colorbar
%    xlabel('\mu*F/0.001');
%    ylabel('F / 0.001');
%    pic_output_color('k1_map');

%    figure(2);
%    set(gca, 'fontsize',font_size);
%    plot(s_prps/0.001, 1000./ISI_a_b(ceil(end/2),:),'-+','linewidth',line_width);
%    xlabel('\mu*F/0.001');
%    ylabel('firing frequency /Hz');
%    pic_output('aveISI');

    fxs = @log;
    inv_fxs = @exp;
    %fxs = inv_tanpatan(s_prps(s_id_prps(1))/0.001, s_prps(s_id_prps(end))/0.001);
    %inv_fxs = tanpatan(s_prps(s_id_prps(1))/0.001, s_prps(s_id_prps(end))/0.001);
    label_num = 6;
    [xl, x_tick, x_tick_val] = xscale_plot(...
        s_prps(s_id_prps)/0.001, [], fxs, inv_fxs, label_num);

    figure(3);  set(gca, 'fontsize',font_size);
    hd = plot(xl, ISI_a_b([1,ceil(size(ISI_a_b,1)/2),end],:),':', xl, 1000./ISI_a_b([1,ceil(size(ISI_a_b,1)/2),end],:),'-+');
    set(hd, 'linewidth',line_width);
    hd=legend(['F=',num2str(s_ps(1)),' ISI'],...
              ['F=',num2str(s_ps(ceil(size(s_ps,2)/2))),' ISI'],...
              ['F=',num2str(s_ps(end)),' ISI'],...
              ['F=',num2str(s_ps(1)),' freq'],...
              ['F=',num2str(s_ps(ceil(size(s_ps,2)/2))),' freq'],...
              ['F=',num2str(s_ps(end)),' freq'],...
              'location','northwest');
    set(hd, 'fontsize',font_size-3);
    if isempty(strfind(signature,'expIF'))
        axis([1.6 3.7 0 200]);
    else
        sa=axis();  sa([3,4])=[0,200];  axis(sa);
    end
    ylabel('ave ISI/ms | Freq/Hz');
    xlabel('\mu*F/0.001');
    set(gca,'xtick', x_tick);
    set(gca,'xticklabel',x_tick_val);
    pic_output_color('aveISI_Freq_xlog');

    figure(5);  set(gca, 'fontsize',font_size);
    [xx2,yy2] = meshgrid(xl, s_ps(s_id_ps)/0.001);
    pcolor(xx2,yy2,k1);
    if mode_st
        caxis([0,100]);
    else
        caxis([0,15]);
    end
    if mode_eif
        ;
    else
        axis([1.6 3.7 0 s_ps(end)/0.001]);
    end
    %shading('flat');
    shading('interp');
    hd = colorbar();
    set(hd, 'fontsize',font_size);
%    title('minGC1/maxGC0 map');
    ylabel('F / 0.001');
    xlabel('\mu*F/0.001');
    set(gca,'xtick', x_tick);
    set(gca,'xticklabel',x_tick_val);
    pic_output_color('k1_map_xlog');
    pic_data_save('k1_map_xlog',xx2,yy2,k1);

    figure(6);  set(gca, 'fontsize',font_size);
    pcolor(xx2,yy2,s_bic_od);
    if mode_st
        caxis([0,100]);
    else
        caxis([0,100]);
    end
    %shading('flat');
    shading('interp');
    hd = colorbar();
    set(hd, 'fontsize',font_size);
    if mode_eif
        ;
    else
        axis([1.6 3.7 0 s_ps(end)/0.001]);
    end
%    title('bic map');
    ylabel('F / 0.001');
    xlabel('\mu*F/0.001');
    set(gca,'xtick', x_tick);
    set(gca,'xticklabel',x_tick_val);
    pic_output_color('bic_map_xlog');

    if mode_st
        caxis_gc_low = 0.1;
        caxis_gc_high= 5.0;
    else
        caxis_gc_low = 0.5;
        caxis_gc_high= 5.0;
    end
    for ii=1:p
    for jj=1:p
        if ii==jj
        continue;
        end
        figure(7);  set(gca, 'fontsize',font_size);
        pcolor(xx2, yy2, s_zero_GC(:,:,ii,jj)/0.001);
        if neu_network(ii,jj)==0
            caxis([0, caxis_gc_low]);
        else
            caxis([0, caxis_gc_high]);
        end
        %shading('flat');
        shading('interp');
        hd = colorbar();
        set(hd, 'fontsize',font_size);
%        title(sprintf('GC%d map %d->%d', neu_network(ii,jj), jj, ii));
        ylabel('F/0.001');
        xlabel('\mu*F/0.001');
        set(gca,'xtick', x_tick);
        set(gca,'xticklabel',x_tick_val);
        pic_name = sprintf('GC%d_map_%d_to_%d_xlog', neu_network(ii,jj), jj, ii);
        pic_output_color(pic_name);
    end
    end

    figure(8);  set(gca, 'fontsize',font_size);
    pcolor(xx2,yy2,wrong_num);
    shading('flat');
    hd = colorbar();
    set(hd, 'fontsize',font_size);
%    if isempty(strfind(signature,'expIF'))
%        axis([1.6 3.7 0 s_ps(end)/0.001]);
%    else
%        ;
%    end
    title(['wrong number (p-value=',num2str(p_value),')']);
    ylabel('F / 0.001');
    xlabel('\mu*F/0.001');
    set(gca,'xtick', x_tick);
    set(gca,'xticklabel',x_tick_val);
    pic_output_color('wrongNum');

end  % stv
end  % scee
end  % simu_time
end  % net

end % id_signature

toc();
