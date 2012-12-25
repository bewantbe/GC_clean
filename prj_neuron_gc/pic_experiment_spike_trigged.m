% Show all spike trigger graph (of 2 neuron)
% only for 2-var
% called by pic_experiment_all

% input:
% X, ras, stv
% linewidth, fontsize, pic_output

varsrd = var(srd,1,2);
varrd  = var(rd,1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% spike-trigged graph
for s_pair=0:3
    p_driving = rem(s_pair,2)+1;
    p_passive = rem(floor(s_pair/2),2)+1;
    st_cond = ['_',char(p_passive+119),'O',char(p_driving+119)];
    ana_len = 200;
    switch p_passive
      case 1
        st_passive = '\epsilon';
      case 2
        st_passive = '\eta';
    end

    %%%%%%%%% volt
    [tg_V, t_rel] = spikeTriggerAve(p_driving, p_passive, ras, X, ana_len, stv);
    figure(4); set(gca, 'fontsize',fontsize);
    plot(0, tg_V(t_rel==0),'rx', t_rel, tg_V, '-+','linewidth',linewidth);
%    axis([min(t_rel), max(t_rel)]);
    xlabel('t_{rel}/ms');
    ylabel(['voltage_{',num2str(p_passive),'|',num2str(p_driving),'}']);
    pic_output(['V', st_cond]);
    pic_data_save(['V', st_cond], 0, tg_V(t_rel==0), t_rel, tg_V);

    %%%%%%%%%% epsilon^*
    [tg_srd, t_rel] = spikeTriggerAve(p_driving, p_passive, ras, srd, ana_len, stv);
    figure(5); set(gca, 'fontsize',fontsize);
    plot(0, tg_srd(t_rel==0),'rx', t_rel, tg_srd, '-+','linewidth',linewidth);
%    axis([min(t_rel), max(t_rel)]);
    xlabel('t_{rel}/ms');
    ylabel([st_passive,'^{*}_{',num2str(p_passive),'|',num2str(p_driving),'}']);
    pic_output(['srd', st_cond]);
    pic_data_save(['srd', st_cond], 0, tg_srd(t_rel==0), t_rel, tg_srd);

    %%%%%%%%%% epsilon
    [tg_rd, t_rel] = spikeTriggerAve(p_driving, p_passive, ras,rd, ana_len, stv);
    figure(6); set(gca, 'fontsize',fontsize);
    plot(0, tg_rd(t_rel==0),'rx', t_rel, tg_rd, '-+','linewidth',linewidth);
%    axis([min(t_rel), max(t_rel)]);
    xlabel('t_{rel}/ms');
    ylabel([st_passive,'_{',num2str(p_passive),'|',num2str(p_driving),'}']);
    pic_output(['rd', st_cond]);
    pic_data_save(['rd', st_cond], 0, tg_rd(t_rel==0), t_rel, tg_rd);

    %%%%%%%%%% epsilon^* - epsilon
    % is this necessary?

    %%%%%%%%%% (epsilon^* - epsilon)^2
    if p_driving ~= p_passive
        tg_gc = (tg_srd - tg_rd).^2 / varrd(p_passive);
        figure(7); set(gca, 'fontsize',fontsize);
        plot(t_rel, tg_gc, 'linewidth',linewidth);
%        axis([min(t_rel), max(t_rel)]);
        xlabel('t_{rel}/ms');
        ylabel(['approx GC_{',num2str(p_passive),'|',num2str(p_driving),'}']);
        pic_output(['tg_gc', st_cond]);
        pic_data_save(['tg_gc', st_cond], t_rel, tg_gc);
    end

    %%%%%%%%%% fire probability
    acce_index = round(ras(ras(:,1)==p_passive, 2) / stv);
    acce_index(acce_index >= len-ana_len | acce_index <= ana_len) = [];
    yf = zeros(1,len);
    yf(acce_index') = 1;
    [tg_fire, t_rel] = spikeTriggerAve(p_driving, p_passive, ras,yf, ana_len, stv);
    figure(8); set(gca, 'fontsize',fontsize);
    plot(0, tg_fire(t_rel==0),'rx', t_rel, tg_fire, '-+','linewidth',linewidth);
%    axis([min(t_rel), max(t_rel)]);
    xlabel('t_{rel}/ms');
    ylabel(['fire probability_{',num2str(p_passive),'|',num2str(p_driving),'}']);
    pic_output(['fire', st_cond]);
    pic_data_save(['fire', st_cond], 0, tg_fire(t_rel==0), t_rel, tg_fire);

end  % x|y
