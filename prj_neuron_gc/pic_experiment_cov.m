% show all cov or cor plots
% only for 2-var
% input  variables: cal_mode, X, srd,  pic_output(), fontsize, linewidth
% output variables: (none)

if ~exist('cal_mode','var')
  cal_mode = 'cov';           % 'cov' or 'cor', calculate corralation or covariance
end
od_show = 40;                           % maximum order to show

%%%%%%%%%%%%%%%%%%        show cov/cor (x, y)       %%%%%%%%%%%%%%%%%%
if strcmp(cal_mode,'cov')
  ssrdR = getcovpd([X(1,:); X(2,:)], od_show);
else
  u0 = X(1,:)*X(1,:)'/len;                        % with bias although
  v0 = X(2,:)*X(2,:)'/len;
  ssrdR = getcovpd([X(1,:)/sqrt(u0); X(2,:)/sqrt(v0)], od_show);
end

%%%%% cov/cor (x(t), x(t-k))
figure(1); set(gca, 'fontsize',fontsize);
plot(1:od_show, ssrdR(1,3:2:end),'o','linewidth',linewidth);
xlabel('k');
ylabel([cal_mode,' (X_{t}, X_{t-k})']);
pic_output([cal_mode,'_x_k']);
pic_data_save([cal_mode,'_x_k'],1:od_show, ssrdR(1,3:2:end));

%%%%% cov/cor (y(t), y(t-k))
figure(2); set(gca, 'fontsize',fontsize);
plot(1:od_show, ssrdR(2,4:2:end),'o','linewidth',linewidth);
xlabel('k');
ylabel([cal_mode,' (Y_{t}, Y_{t-k})']);
pic_output([cal_mode,'_y_k']);
pic_data_save([cal_mode,'_y_k'],1:od_show, ssrdR(2,4:2:end));

%%%%% cov/cor (x(t), y(t-k))
figure(3); set(gca, 'fontsize',fontsize);
plot(0,ssrdR(1,2),'*',-od_show:od_show, [fliplr(ssrdR(2,3:2:end)), ssrdR(1,2:2:end)],'o','linewidth',linewidth);
xlabel('k');
ylabel([cal_mode,' (X_{t}, Y_{t-k})']);
pic_output([cal_mode,'_xy_k']);
pic_data_save([cal_mode,'_xy_k'],0,ssrdR(1,2),'*',-od_show:od_show, [fliplr(ssrdR(2,3:2:end)), ssrdR(1,2:2:end)]);

%%%%%%%%%%%%%%%%%%   show cov/cor (\epsilon,\eta)   %%%%%%%%%%%%%%%%%%
if strcmp(cal_mode,'cov')
  ssrdR = getcovpd([srd(1,:); srd(2,:)], od_show);
else
  u0 = srd(1,:)*srd(1,:)'/len;                   % with bias although
  v0 = srd(2,:)*srd(2,:)'/len;
  ssrdR = getcovpd([srd(1,:)/sqrt(u0); srd(2,:)/sqrt(v0)], od_show);
end

%%%%% cov/cor (epsilon^*(t), epsilon^*(t-k))
figure(1); set(gca, 'fontsize',fontsize);
plot(1:od_show, ssrdR(1,3:2:end),'o','linewidth',linewidth);
xlabel('k');
ylabel([cal_mode,' (\epsilon^{*}_{t}, \epsilon^{*}_{t-k})']);
pic_output([cal_mode,'_srdx_k']);
pic_data_save([cal_mode,'_srdx_k'],1:od_show, ssrdR(1,3:2:end));

%%%%% cov/cor (eta^*(t), eta^*(t-k))
figure(2); set(gca, 'fontsize',fontsize);
plot(1:od_show, ssrdR(2,4:2:end),'o','linewidth',linewidth);
xlabel('k');
ylabel([cal_mode,' (\eta^{*}_{t}, \eta^{*}_{t-k})']);
pic_output([cal_mode,'_srdy_k']);
pic_data_save([cal_mode,'_srdy_k'],1:od_show, ssrdR(2,4:2:end));

%%%%% cov/cor (epsilon^*(t), eta^*(t-k))
figure(3); set(gca, 'fontsize',fontsize);
plot(0,ssrdR(1,2),'*', -od_show:od_show, [fliplr(ssrdR(2,3:2:end)), ssrdR(1,2:2:end)],'-o','linewidth',linewidth);
xlabel('k');
ylabel([cal_mode,' (\epsilon^{*}_{t}, \eta^{*}_{t-k})']);
pic_output([cal_mode,'_srdxy_k']);
pic_data_save([cal_mode,'_srdxy_k'], 0,ssrdR(1,2),'*', -od_show:od_show, [fliplr(ssrdR(2,3:2:end)), ssrdR(1,2:2:end)]);

%%% for comparison (to spike trigger)
od_show = 75;
if strcmp(cal_mode,'cov')
  ssrdR = getcovpd([srd(1,:); srd(2,:)], od_show);
else
  ssrdR = getcovpd([srd(1,:)/sqrt(u0); srd(2,:)/sqrt(v0)], od_show);
end
%save('-v7','ssrdR.mat','ssrdR');
%load('ssrdR.mat');
%%%% cov/cor (epsilon^*(t), eta^*(t-k))
figure(3); set(gca, 'fontsize',fontsize);
rlr     = -fliplr([fliplr(ssrdR(2,3:2:end)), ssrdR(1,2:2:end)]);
rg      = -od_show:od_show;
rg_show = round(od_show-od_show/3+1):2*od_show+1; % -25:75
hd=plot(rg(rg_show), rlr(rg_show),'k-o', 0,rlr(rg==0),'b*');
set(hd,'linewidth',linewidth);
xlabel('-k');
ylabel(['-',cal_mode,' (\epsilon^{*}_{t}, \eta^{*}_{t-k})']);
pic_output_color([cal_mode,'_-srdxy_-k']);
pic_data_save([cal_mode,'_-srdxy_-k'], rg(rg_show), rlr(rg_show), 0,rlr(rg==0));

