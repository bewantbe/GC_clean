% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat
% constant data length
t0=tic();
tocs = @(st) fprintf('%s: t = %6.3fs\n', st, toc());

% _ST _expIF
signature = 'data_scan_stv/net_2_2_IF_sc2_l1e7_w21';  % to distinguish different parallel program instances (also dir)

if isempty(strfind(upper(signature),upper('expIF')))
  %s_prps_default = logspace(log10(4.9e-3), log10(4.7e-2), 30);
  dt_std = 1.0/32;
  model_mode = 'IF';
else
  %s_prps_default = logspace(log10(1.5e-3), log10(1.2e-1), 30);
  dt_std = 0.004;
  model_mode = 'EIF';
end
b_spike_train = isempty(strfind(upper(signature),upper('ST')))==1;  % using Spike Train?
%s_ps_default = [0.001:0.001:0.029];  %29

% scan value sets
s_net  = {'net_2_2'};
s_time = [1e7];
s_scee = [0.02];
s_prps = [0.012];
s_ps   = [0.012];
s_stv  = [0.5:0.5:20];
maxod  = 99;
s_od   = 1:maxod;
hist_div = 0:0.5:400;         % ISI
T_segment = 1000;             % in ms
stv0   = 0.125;               % fine sample rate

save('-v7', [signature, '_info.mat'], 's_net', 's_time', 's_scee', 's_prps', 's_ps', 's_stv', 's_od', 'hist_div', 'maxod', 'T_segment','stv0');

data_path = ['data/', signature, '_'];

for net_id = 1:length(s_net)
 netstr = s_net{net_id};
 matname = ['network/', netstr, '.txt'];
 neu_network = load('-ascii', matname);
 p = size(neu_network, 1);
for simu_time = s_time
for scee = s_scee
 prps_ps_stv_oGC = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_Sgc = zeros(p, p, length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_oDe = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_Sde = zeros(p, length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_R   = zeros(p, p*(maxod+1), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_S   = cell(length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_fqs = cell(length(s_prps), length(s_ps), length(s_stv));
 prps_ps_aveISI  = zeros(p, length(s_prps), length(s_ps));
 prps_ps_ISI_dis = zeros(p, length(hist_div), length(s_prps), length(s_ps));
for id_prps = 1:length(s_prps)
 prps = s_prps(id_prps);
for id_ps = 1:length(s_ps)
 ps = s_ps(id_ps);
 pr = prps / ps;
for id_stv = 1:length(s_stv)
    stv = s_stv(id_stv);
    disp(sprintf('ps %f, pr %f, prps %f, stv %f',ps,pr,prps,stv));
    if exist('OCTAVE_VERSION','builtin')
        fflush(stdout);
    end
    if id_stv == 1     % actuall calculation
      extpara = '--RC-filter';
      [oX, aveISI, ras] = gendata_neu(netstr, scee, pr, ps, simu_time, stv0, extpara);
      [p, len] = size(oX);
      mlen = T_segment/stv0;
      % read volt data, binary format
      if b_spike_train
          % get spike train
          len = round(simu_time/s_stv(1));
          oX = zeros(p,len);
          for neuron_id = 1:p
              st = SpikeTrain(ras, len, neuron_id, 1, stv);
              oX(neuron_id,:) = st;
          end
      end
    end

    slen = round(T_segment/stv);
    [X1,X2,T1,T2] = SampleNonUnif(oX, mlen, slen, '1', simu_time/T_segment*stv/s_stv(1));
    fqs = ifftshift(((0:slen-1)-floor(slen/2))/slen);

    tic();
    S_x2 = mX2S_nuft(X1, T1);
    S1 = Makeup4SpectrumFact(S_x2);
    [gc, de11, de22] = getGCSapp(S1);
    tocs('time spend in GC app');

    [oGC, oDe, R] = AnalyseSeries(oX(:, 1:(stv/stv0):end), s_od);

    prps_ps_stv_oGC(:,:,:, id_prps, id_ps, id_stv) = oGC;
    prps_ps_stv_Sgc(:, :, id_prps, id_ps, id_stv) = gc;
    prps_ps_stv_oDe(:,:,:, id_prps, id_ps, id_stv) = oDe;
    prps_ps_stv_Sde(:, id_prps, id_ps, id_stv) = [de11, de22];
    prps_ps_stv_R  (:,:, id_prps, id_ps, id_stv) = R;
    prps_ps_stv_S  {id_prps, id_ps, id_stv} = S1;
    prps_ps_stv_fqs{id_prps, id_ps, id_stv} = fqs;
end  % stv
    prps_ps_aveISI(:, id_prps, id_ps) = aveISI;
    % get ISI distribution
    ISI_dis = zeros(p, length(hist_div));
    for kk = 1:p
        tm = ras(ras(:,1)==kk, 2);
        tm = tm(2:end) - tm(1:end-1);
        ISI_dis(kk,:) = hist(tm, hist_div);
    end
    prps_ps_ISI_dis(:,:,id_prps, id_ps) = ISI_dis;
end  % ps
end  % prps
 datamatname = sprintf('%s_%s_sc=%g_t=%.3e.mat', signature, netstr, scee, simu_time);
 save('-v7', datamatname, 'prps_ps_stv_oGC', 'prps_ps_stv_Sgc', 'prps_ps_stv_oDe', 'prps_ps_stv_Sde', 'prps_ps_stv_fqs', 'prps_ps_stv_R', 'prps_ps_stv_S', 'prps_ps_aveISI', 'prps_ps_ISI_dis');
 % save length ?
end  % scee
end  % simu_time
end  % net

fprintf('Elapsed time is %6.3f\n', (double(tic()) - double(t0))*1e-6 );
% vim: set ts=4 sw=4 ss=4
