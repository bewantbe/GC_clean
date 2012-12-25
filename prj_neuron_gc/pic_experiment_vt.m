% show voltage trace (time course)
% input  variables: rg, stv, X, srd, rd, fontsize, linewidth 
% output variables: (none)

% 'plump'
% id_ed=11500; rg = id_ed-199:id_ed; pic_experiment_vt
old_fontsize = fontsize;
%fontsize = 18;
fontsize = 24;

% plot range
%if ~exist('rg','var')
%  disp('warning: rg not defined! (pic_experiment_vt)');
%  id_ed = 11500;
%  rg    = id_ed-199:id_ed;
%end
id_ed = 11500;
rg    = id_ed-199:id_ed;
if strcmpi(cal_model,'middle2')
  id_ed = 11500;
  rg    = id_ed-199:id_ed;
end
if ~isempty(strcmpi(cal_model,'_sync'))
  id_ed = 11500;
  rg    = id_ed-1000:id_ed-1;
end
if p==2
  % voltage
  figure(6);  set(gca, 'fontsize',fontsize);
  hd=plot((rg-rg(1))*stv, X(1,rg),(rg-rg(1))*stv, X(2,rg),'--');
  set(hd, 'linewidth',linewidth);
  ylabel('voltage');
  xlabel('t /ms');
  hd=legend('X','Y');
%  hd=legend('right');
  set(hd, 'fontsize',fontsize);
  pic_output_color('volt_t');
  pic_data_save('volt_t',(rg-rg(1))*stv, X(1,rg),(rg-rg(1))*stv, X(2,rg));

  % show residual
  figure(7);  set(gca, 'fontsize',fontsize);
  hd=plot((rg-rg(1))*stv, srd(1,rg),(rg-rg(1))*stv, srd(2,rg),'--','linewidth',linewidth);
  set(hd, 'linewidth',linewidth);
  ylabel('auto-regression residual');
  xlabel('t /ms');
  hd=legend('X','Y');
%  hd=legend('right');
  set(hd, 'fontsize',fontsize);
  pic_output_color('srd_t');
  pic_data_save('srd_t',(rg-rg(1))*stv, srd(1,rg),(rg-rg(1))*stv, srd(2,rg));

  figure(8);  set(gca, 'fontsize',fontsize);
  hd=plot((rg-rg(1))*stv,  rd(1,rg),(rg-rg(1))*stv,  rd(2,rg),'--','linewidth',linewidth);
  set(hd, 'linewidth',linewidth);
  ylabel('joint-regression residual');
  xlabel('t /ms');
  hd=legend('X','Y');
%  hd=legend('right');
  set(hd, 'fontsize',fontsize);
  pic_output_color('rd_t');
  pic_data_save('rd_t',(rg-rg(1))*stv,  rd(1,rg),(rg-rg(1))*stv,  rd(2,rg));

  figure(9);  set(gca, 'fontsize',fontsize);
  tmp_gct = srd(:,rg) - rd(:,rg);
  hd=plot((rg-rg(1))*stv, tmp_gct(1,:),(rg-rg(1))*stv, tmp_gct(2,:),'--','linewidth',linewidth);
  xlabel('t /ms');
  hd=legend('\epsilon^*_t - \epsilon_t', '\eta^*_t - \eta_t');
%  hd=legend('right');
  set(hd, 'fontsize',fontsize);
  pic_output_color('[srd-rd]_t');
  pic_data_save('[srd-rd]_t',(rg-rg(1))*stv, tmp_gct(1,:),(rg-rg(1))*stv, tmp_gct(2,:));

  tmp_gct = (srd(:,rg) - rd(:,rg)).^2;
  figure(10); set(gca, 'fontsize',fontsize);
  hd=plot((rg-rg(1))*stv, tmp_gct(1,:), (rg-rg(1))*stv, tmp_gct(2,:),'--','linewidth',linewidth);
  set(hd, 'linewidth',linewidth);
  ylabel('GC time expansion');
  xlabel('t /ms');
  hd=legend('(\epsilon^*_t - \epsilon_t)^2, Y->X', '(\eta^*_t - \eta_t)^2, X->Y');
%  hd=legend('right');
  set(hd, 'fontsize',fontsize);
  pic_output_color('[srd-rd]^2_t');
  pic_data_save('[srd-rd]^2_t',(rg-rg(1))*stv, tmp_gct(1,:), (rg-rg(1))*stv, tmp_gct(2,:));
else % p > 2
  id_legend = cell(1,p);
  for k=1:p
    id_legend{k} = num2str(k);
  end
  % voltage
  figure(6);  set(gca, 'fontsize',fontsize);
  plot((rg-rg(1))*stv, X(:,rg),'linewidth',linewidth);
  ylabel('voltage');
  xlabel('t /ms');
  hd=legend(id_legend);  set(hd, 'fontsize',fontsize);
  pic_output_color('volt_t');
  pic_data_save('volt_t',(rg-rg(1))*stv, X(:,rg));

  % show residual
  figure(7);  set(gca, 'fontsize',fontsize);
  plot((rg-rg(1))*stv, srd(:,rg),'linewidth',linewidth);
  ylabel('auto-regression residual');
  xlabel('t /ms');
  hd=legend(id_legend);  set(hd, 'fontsize',fontsize);
  pic_output_color('srd_t');
  pic_data_save('srd_t',(rg-rg(1))*stv, srd(:,rg));

  figure(8);  set(gca, 'fontsize',fontsize);
  plot((rg-rg(1))*stv, rd(:,rg),'linewidth',linewidth);
  ylabel('joint-regression residual');
  xlabel('t /ms');
  hd=legend(id_legend);  set(hd, 'fontsize',fontsize);
  pic_output_color('rd_t');
  pic_data_save('rd_t',(rg-rg(1))*stv, rd(:,rg));

  figure(9);  set(gca, 'fontsize',fontsize);
  plot((rg-rg(1))*stv, srd(:,rg) - rd(:,rg),'linewidth',linewidth);
  ylabel('\epsilon^*_t - \epsilon_t');
  xlabel('t /ms');
  hd=legend(id_legend);  set(hd, 'fontsize',fontsize);
  pic_output_color('[srd-rd]_t');
  pic_data_save('[srd-rd]_t',(rg-rg(1))*stv, srd(:,rg) - rd(:,rg));

  figure(10); set(gca, 'fontsize',fontsize);
  plot((rg-rg(1))*stv, (srd(:,rg) - rd(:,rg)).^2,'linewidth',linewidth);
  ylabel('(\epsilon^*_t - \epsilon_t)^2');
  xlabel('t /ms');
  hd=legend(id_legend);  set(hd, 'fontsize',fontsize);
  pic_output_color('[srd-rd]^2_t');
  pic_data_save('[srd-rd]^2_t',(rg-rg(1))*stv, (srd(:,rg) - rd(:,rg)).^2);
end
fontsize = old_fontsize;
