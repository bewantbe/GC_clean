% Show all analysis graph

% whether we are in Octave or Matlab
is_octave = exist('OCTAVE_VERSION','builtin') ~= 0;

if is_octave
  fontsize = 28;
  linewidth = 5;
else
  fontsize = 24;
  linewidth = 3;
end

% fixed parameters
netstr = 'net_2_2';
scee = 0.01;
%pr = 1.7;  % pr = 0.24, 0.38, 0.5, 0.8, 1.7; ps=0.02,scee=0.01,simu_time=1e7 or 5e7 or 2e7
ps = 0.02;
simu_time = 5e5;                             % 1e7 or 5e7
stv = 0.5;
%mode_st = true;
mode_st = false;
mode_eif = false;

% default values
%pic_prefix0 = ['pic_spiketrigged_matlab_5e7/',netstr,'_para01hp'];
%s_cal_model = {'plump','inflammable','middle','middle2'};  % cost 8059.5sec using MATLAB R2011b in xhp; but only 2025.7sec in octave with non-optimised ATLAS
% cost 530.99 sec at hi7
%pic_prefix0 = ['pic_spiketrigged/',netstr];
%s_cal_model = {'test1_sync', 'test2_sync', 'test3_sync'};

%% overwrite default values
pic_prefix0 = ['pic_tmp/',netstr,'_para01hp'];
s_cal_model = {'plump','inflammable','middle','middle2'};

if mode_eif
  pic_prefix0 = [pic_prefix0, '_EIF'];
end

if mode_st
  pic_prefix0 = [pic_prefix0, '_ST'];
end

tic;
%%%%% big for loop
for id_cal_model=2:length(s_cal_model)
cal_model = s_cal_model{id_cal_model};
pic_prefix = [pic_prefix0, '_', cal_model, '_'];
%pic_output = @(st)print('-dpng',[pic_prefix, st, '.png'],'-r100');    % output function
%pic_output_color = pic_output;                                        % for color output
pic_output = @(st)print('-deps',[pic_prefix, st, '.eps']);
pic_output_color = @(st)print('-depsc2',[pic_prefix, st, '.eps']);

% for saving data of each picture
%   vx1=varargin{1}, vy1=varargin{2}, vx2=varargin{2} ...
%pic_data_save = @(st, vx, vy)save('-v7', [pic_prefix, st, '.mat'], 'vx', 'vy');
pic_data_save = @(st, varargin)save('-v7', [pic_prefix, st, '.mat'], 'varargin');

clear('X','srd','rd');                       % save memory
%%%%%%%%%%%% get data
if mode_eif
  ext_cmd = 'ExpIF -dt 0.004 --RC-filter';
else
  ext_cmd = '-dt 0.0078125 --RC-filter';
end
if strcmp(cal_model, 'plump') == 1
  pr = 0.24;  ps = 0.02;
end
if strcmp(cal_model, 'inflammable') == 1
  pr = 1.7;  ps = 0.02;
end
if strcmp(cal_model, 'middle') == 1
  pr = 0.38;  ps = 0.02;
end
if strcmp(cal_model, 'middle2') == 1
  pr = 0.5;  ps = 0.02;
end
if strcmp(cal_model, 'test1_sync') == 1
  pr = 7;  ps = 0.001;  scee = 0.01;
end
if strcmp(cal_model, 'test2_sync') == 1
  pr = 7;  ps = 0.001;  scee = 0.01;
  ext_cmd = 'new -dt 0.0078125 --RC-filter --pr-mul 1.031 1';  % fr=19.69 19.69
end
if strcmp(cal_model, 'test3_sync') == 1
  pr = 7;  ps = 0.001;  scee = 0.02;
  ext_cmd = 'new -dt 0.0078125 --RC-filter --pr-mul 1.07142857 1';  % fr=23.11
end
[X, ISI, ras] = gendata_neu(netstr, scee, pr, ps, simu_time, stv, ext_cmd);
[p, len] = size(X);
if mode_st
    for neuron_id = 1:p
        st = SpikeTrain(ras, len, neuron_id, 1, stv);
        X(neuron_id,:) = st;
    end
end

disp('ISI:');
disp(ISI);
fprintf('net:%s, sc:%.3f, pr:%.2f, ps:%.4f, time:%.2e, stv:%.2f, len:%.2e\n',...
 netstr, scee, pr, ps, simu_time, stv, len);

%%%%%%%%%%%% analyse
if mode_st
  od_max = 199;
else
  od_max = 100;
end
if mode_eif
  od_max = 499;
end
GC_regression;

fprintf('bo:%2d, ao:%2d\n', bic_od, aic_od);
disp('var (\Sigma) at bic_od:');
disp(oDe(:,:,bic_od));
format 'bank'
  disp('GC(od=20)*1e4:');
  disp(1e4*oGC(:,:,20));
  disp(['GC(od=', num2str(bic_od), ')*1e4:']);
  disp(1e4*oGC(:,:,bic_od));
  disp('GC(od=zero)*1e4:');
  disp(1e4*zero_GC);
format 'short'
  disp('non-zero test:');
  disp(gc_prob_nonzero(oGC(:,:,bic_od),bic_od,len));
disp('var ratio (at bic_od): var(x)/var(\epsilon), var(y)/var(\eta)');   % for 2-var only
disp([num2str(R(1,1)/oDe(1,1,bic_od)), ' ', num2str(R(2,2)/oDe(2,2,bic_od))]);

%%%%%%%%%% coef of regression of epsilon^*
od_show = 20;
[rdgc, rdde, rdA] = pos_nGrangerT2(srd, bic_od);
od = min([od_show, bic_od]);
figure(9);  set(gca, 'fontsize',fontsize);
plot(1:od, rdA(2,1:2:2*od),'*','linewidth',linewidth);
xlabel('j');
ylabel('c_j');
pic_output('c_j');
pic_data_save('c_j',1:od,rdA(2,1:2:2*od));

figure(10); set(gca, 'fontsize',fontsize);
plot(1:od, rdA(1,2:2:2*od),'*','linewidth',linewidth);
xlabel('j');
ylabel('b_j');
pic_output('b_j');
pic_data_save('b_j',1:od,rdA(1,2:2:2*od));

pic_experiment_vt;

pic_experiment_pdf;

cal_mode = 'cov';
pic_experiment_cov;

cal_mode = 'cor';
pic_experiment_cov;

pic_experiment_spike_trigged;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine
% <pic name='coef: c_{j}'>
% <pic name='cross-correlation of \eps and \eta k=1~20'>
od_show = 30;
od_show = min([od_show, bic_od]);
[rdgc, rdde, rdA] = pos_nGrangerT2(srd, bic_od);
ssrdR = getcovpd([srd(1,:); srd(2,:)], od_show);

figure(10);  set(gca, 'fontsize',fontsize);
%plot(1:od_show, rdA(2,1:2:2*od_show),'+','markersize',fontsize*0.5, 1:od_show,-ssrdR(2,3:2:end)/var(rd(1,:)),'o','markersize',fontsize*0.5);
plot(1:od_show, rdA(2,1:2:2*od_show),'+', 1:od_show,-ssrdR(2,3:2:end)/var(rd(1,:)),'o','markersize',fontsize*0.5);
xlabel('j');
hd=legend('c_j','-cov(\epsilon^{*}_{t},\eta^{*}_{t+k})/var(\epsilon^{*}_{t})','location','southeast');
%hd=legend('right');
set(hd,'fontsize',fontsize-4);
pic_output('c_j_and_cov_srdyx_j');
pic_data_save('c_j_and_cov_srdyx_j',1:od_show, rdA(2,1:2:2*od_show), 1:od_show,-ssrdR(2,3:2:end)/var(rd(1,:)));

figure(11);  set(gca, 'fontsize',fontsize);
plot(1:od_show, rdA(2,1:2:2*od_show)+ssrdR(2,3:2:end)/var(rd(1,:)),'+','markersize',fontsize*0.5);
xlabel('j');
ylabel('c_j+cov(\epsilon^{*}_{t},\eta^{*}_{t+k})/var(\epsilon^{*}_{t})');
pic_output('c_j_and_cov_srdyx_j_diff');
pic_data_save('c_j_and_cov_srdyx_j_diff',1:od_show, rdA(2,1:2:2*od_show)+ssrdR(2,3:2:end)/var(rd(1,:)));

end  % big parameter for loop

toc

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine
% <pic name='coef: c_{j}'>
% <pic name='cross-correlation of \eps and \eta k=1~20'>
od_show = 30;
od_show = min([od_show, bic_od]);
[rdgc, rdde, rdA] = pos_nGrangerT2(srd, bic_od);
ssrdR = getcovpd([srd(1,:); srd(2,:)], od_show);

c_max = max(-rdA(2,1:2:2*od_show));
r_max = max(ssrdR(2,3:2:end));

figure(10);  set(gca, 'fontsize',fontsize);
%plotyy(1:od_show,-rdA(2,1:2:2*od_show), 1:od_show,ssrdR(2,3:2:end));
plot(1:od_show,-rdA(2,1:2:2*od_show)/c_max,'+','markersize',fontsize*0.5, 1:od_show,ssrdR(2,3:2:end)/r_max,'o','markersize',fontsize*0.5);
xlabel('j');
ylabel('normalized scale');
hd=legend('-c_j','cor(\epsilon^{*}_{t}, \eta^{*}_{t+k})');
set(hd,'fontsize',fontsize-4);
pic_output('c_j_and_cov_srdyx_j');

figure(11);  set(gca, 'fontsize',fontsize);
plot(1:od_show,-rdA(2,1:2:2*od_show)/c_max-ssrdR(2,3:2:end)/r_max,'+','markersize',fontsize*0.5);
xlabel('j');
ylabel('normalized -c_j-cor(\epsilon^{*}_{t}, \eta^{*}_{t+k})');
pic_output('c_j_and_cov_srdyx_j_diff');
%}

