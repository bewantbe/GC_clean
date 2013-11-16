% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat
% constant data length
tic();

% _ST _expIF
signature = 'data_scan/template_w01';     % to distinguish different parallel program instances (also dir)

if isempty(strfind(upper(signature),upper('expIF')))
  s_prps_default = logspace(log10(4.9e-3), log10(4.7e-2), 30);
  dt_std = 1.0/32;
  static_param = ['../raster_tuning       -ng -q --bin-save --RC-filter -dt ', num2str(dt_std,'%.16g')];
else
  s_prps_default = logspace(log10(1.5e-3), log10(1.2e-1), 30);
  dt_std = 0.004;
  static_param = ['../raster_tuning_expIF -ng -q --bin-save --RC-filter -dt ', num2str(dt_std,'%.16g')];
end
b_spike_train = isempty(strfind(upper(signature),upper('ST')))==0;  % using Spike Train?
s_ps_default = [0.001:0.001:0.029];  %29

% scan value sets
s_net  = {'net_3_06'};
s_time = [1e5];
s_scee = [0.02];
s_prps = s_prps_default;  % 30
s_ps   = s_ps_default;    % 29
s_stv  = [0.5];
maxod  = 499;
s_od   = 1:maxod;
hist_div = 0:0.5:400;          % ISI

save('-v7', [signature, '_info.mat'], 's_net', 's_time', 's_scee', 's_prps', 's_ps', 's_stv', 's_od', 'hist_div', 'maxod');

data_path = ['data/', signature, '_'];

for net_id = 1:length(s_net)
 netstr = s_net{net_id};
 matname = ['network/', netstr, '.txt'];
 neu_network = load('-ascii', matname);
 p = size(neu_network, 1);
 output_name     = [data_path, 'volt_', netstr, '.dat'];
 output_ISI_name = [data_path, 'ISI_',  netstr, '.txt'];
 output_RAS_name = [data_path, 'RAS_',  netstr, '.txt'];
for simu_time = s_time
for scee = s_scee
 prps_ps_stv_oGC = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_oDe = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_R   = zeros(p, p*(maxod+1), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_aveISI  = zeros(p, length(s_prps), length(s_ps));
 prps_ps_ISI_dis = zeros(p, length(hist_div), length(s_prps), length(s_ps));
for id_prps = 1:length(s_prps)
 prps = s_prps(id_prps);
for id_ps = 1:length(s_ps)
 ps = s_ps(id_ps);
 pr = prps / ps;
 disp(sprintf('ps %f, pr %f, prps %f',ps,pr,prps));
 if exist('octave_config_info','builtin')
   fflush(stdout);
 end
 cmdst_neu = sprintf('-n %d -mat %s -scee %.16e -pr %.16e -ps %.16e', ...
                     p, matname, scee, pr, ps);
for id_stv = 1:length(s_stv)
    stv = s_stv(id_stv);
    if id_stv == 1                 % output the first data, since it's longer
        cmdst_simu= sprintf('-t %.16e -stv %.16e -o "%s" --save-spike-interval %s --save-spike %s', ...
                     simu_time, stv, output_name, output_ISI_name, output_RAS_name);
    else
        cmdst_simu= sprintf('-t %.16e -stv %.16e -o "%s"', ...
                     simu_time, stv, output_name);
    end
    system([static_param,' ',cmdst_neu,' ',cmdst_simu]);

    % read volt data, binary format
    if b_spike_train
        % get spike train
        ras = load('-ascii', output_RAS_name);
        len = round(simu_time/s_stv(1));
        X = zeros(p,len);
        for neuron_id = 1:p
            st = SpikeTrain(ras, len, neuron_id, 1, stv);
            X(neuron_id,:) = st;
        end
    else
        fid = fopen(output_name, 'r');
        X = fread(fid, [p, Inf], 'double');
        fclose(fid);
        [p, len] = size(X);
    end

    [oGC, oDe, R] = AnalyseSeries(X, s_od);

    prps_ps_stv_oGC(:,:,:, id_prps, id_ps, id_stv) = oGC;
    prps_ps_stv_oDe(:,:,:, id_prps, id_ps, id_stv) = oDe;
    prps_ps_stv_R  (:,:, id_prps, id_ps, id_stv) = R;
end  % stv
    prps_ps_aveISI(:, id_prps, id_ps) = load('-ascii', output_ISI_name);
    % get ISI distribution
    ISI_dis = zeros(p, length(hist_div));
    ras = load('-ascii', output_RAS_name);
    for kk = 1:p
        tm = ras(ras(:,1)==kk, 2);
        tm = tm(2:end) - tm(1:end-1);
        ISI_dis(kk,:) = hist(tm, hist_div);
    end
    prps_ps_ISI_dis(:,:,id_prps, id_ps) = ISI_dis;
  toc();
end  % ps
end  % prps
 datamatname = sprintf('%s_%s_sc=%g_t=%.3e.mat', signature, netstr, scee, simu_time);
 save('-v7', datamatname, 'prps_ps_stv_oGC', 'prps_ps_stv_oDe', 'prps_ps_stv_R', 'prps_ps_aveISI', 'prps_ps_ISI_dis');
 % save length ?
end  % scee
end  % simu_time
end  % net

toc();
