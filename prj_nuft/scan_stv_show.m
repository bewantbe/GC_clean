% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat
% constant data length
tic();

% data
%signature = 'data_scan_stv/net_2_2_IF_sc2_l1e7_w3';  % xhp: 31960.818 sec
%signature = 'data_scan_stv/net_2_2_IF_sc2_l1e7_w4';  % 32099 sec
%signature = 'data_scan_stv/net_2_2_IF_sc2_l1e7_w12';
%signature = 'data_scan_stv/net_2_2_IF_sc2_l1e7_w11';  % 32099 sec
%signature = 'data_scan_stv/net_2_2_IF_sc2_l1e7_w21';  % 16644 sec
signature = 'data_scan_stv/net_2_2_IF_sc2_l1e8_w31';  % 1049 sec

if isempty(strfind(upper(signature),upper('expIF')))
  %s_prps_default = logspace(log10(4.9e-3), log10(4.7e-2), 30);
  dt_std = 1.0/32;
  model_mode = 'IF';
  %static_param = ['../raster_tuning       -ng -q --bin-save --RC-filter -dt ', num2str(dt_std,'%.16g')];
else
  %s_prps_default = logspace(log10(1.5e-3), log10(1.2e-1), 30);
  dt_std = 0.004;
  model_mode = 'EIF';
  %static_param = ['../raster_tuning_expIF -ng -q --bin-save --RC-filter -dt ', num2str(dt_std,'%.16g')];
end
b_spike_train = isempty(strfind(upper(signature),upper('ST')))==1;  % using Spike Train?
%s_ps_default = [0.001:0.001:0.029];  %29

% scan value sets
s_net  = {'net_2_2'};
s_time = [2e7];
s_scee = [0.02];
s_prps = [0.012];
s_ps   = [0.012];
s_stv  = [0.5:0.5:1];
maxod  = 99;
s_od   = 1:maxod;
hist_div = 0:0.5:400;         % ISI
T_segment = 64;               % in ms
stv0   = 0.125;                   % fine sample rate

%load 's_net', 's_time', 's_scee', 's_prps', 's_ps', 's_stv', 's_od', 'hist_div', 'maxod', 'T_segment','stv0'
load([signature, '_info.mat']);

data_path = ['data/', signature, '_'];

for net_id = 1:length(s_net)
 netstr = s_net{net_id};
 matname = ['network/', netstr, '.txt'];
 neu_network = load('-ascii', matname);
 p = size(neu_network, 1);
for simu_time = s_time
for scee = s_scee
 datamatname = sprintf('%s_%s_sc=%g_t=%.3e.mat', signature, netstr, scee, simu_time);
 %load 'prps_ps_stv_oGC', 'prps_ps_stv_Sgc', 'prps_ps_stv_oDe', 'prps_ps_stv_Sde', 'prps_ps_stv_R', 'prps_ps_stv_S', 'prps_ps_stv_fqs', 'prps_ps_aveISI', 'prps_ps_ISI_dis'
 load(datamatname);
for id_prps = 1:length(s_prps)
 prps = s_prps(id_prps);
for id_ps = 1:length(s_ps)
 ps = s_ps(id_ps);
 pr = prps / ps;
 s_gc = zeros(p*p-p,length(s_stv));
for id_stv = 1:length(s_stv)
    stv = s_stv(id_stv);

    od = 20;
    %prps_ps_stv_oGC(:,:,:, id_prps, id_ps, id_stv) = oGC;
    %prps_ps_stv_R  (:,:, id_prps, id_ps, id_stv) = R;
    S1 = prps_ps_stv_S  {id_prps, id_ps, id_stv};
    fqs= prps_ps_stv_fqs{id_prps, id_ps, id_stv};

    %gc = prps_ps_stv_oGC(:,:,od, id_prps, id_ps, id_stv);
    %de = prps_ps_stv_oDe(:,:,od, id_prps, id_ps, id_stv);
    %s_Sde(:,id_stv) = de(1:(p+1):(p*p));

    gc = prps_ps_stv_Sgc(:, :, id_prps, id_ps, id_stv);
    s_Sde(:,id_stv) = prps_ps_stv_Sde(:, id_prps, id_ps, id_stv);
    gc(eye(p)==1) = [];
    s_gc_p(:,id_stv) = reshape(gc,[],1);

    %figure(id_stv+4);
    %plot(fftshift(fqs), fftshift(S1(:,1,1)));
end  % stv
    %prps_ps_aveISI(:, id_prps, id_ps) = aveISI;
    %prps_ps_ISI_dis(:,:,id_prps, id_ps) = ISI_dis;

    % save pics
    figure(1);
    plot(s_stv, s_gc_p(:,:),'-o');
    legend('1->2', '2->1');
    figure(2);
    plot(s_stv, s_Sde(:,:));
end  % ps
end  % prps
     %load
end  % scee
end  % simu_time
end  % net

toc();
% vim: ts=4 sw=4 ss=4
