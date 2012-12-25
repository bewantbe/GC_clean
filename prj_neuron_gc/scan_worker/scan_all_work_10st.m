% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat
% constant data length
tic();
addpath ..

signature = 'data_scan_IF/w_net_3_06_sc2_t1e6_199_st';     % to distinguish different parallel program instances (also dir)

% scan value sets
s_net  = {'net_3_06'};
s_time = [1e6];
s_scee = [0.02];
s_prps = [0.005:0.0001:0.00725, 0.0075:0.0005:0.0095, 0.01:0.002:0.019, 0.02:0.005:0.04];  % 38
s_ps   = [0.001:0.0005:0.0055, 0.006:0.001:0.029, 0.03:0.002:0.04];  %40
s_stv  = [0.5];  s_dt   = 1.0/32;
maxod  = 199;
s_od   = 1:maxod;
hist_div = 0:0.5:400;          % ISI

save('-v7', [signature, '_info.mat'], 's_net', 's_time', 's_scee', 's_prps', 's_ps', 's_stv', 's_od', 'hist_div', 'maxod');

static_param = ['../raster_tuning -ng -q --bin-save --RC-filter -dt ', num2str(s_dt(1))];
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
 len = round(simu_time/s_stv(1));
 X = zeros(p,len);
for scee = s_scee
 prps_ps_stv_oGC = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_oDe = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_R   = zeros(p, p*(maxod+1), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_aveISI  = zeros(p, length(s_prps), length(s_ps));
 prps_ps_ISI_dis = zeros(p, length(hist_div), length(s_prps), length(s_ps));
 id_prps = 0;
for prps = s_prps
 id_prps = id_prps + 1;
 id_ps = 0;
for ps = s_ps
 id_ps = id_ps + 1;
 pr = prps / ps;
 disp(sprintf('ps %f, pr %f, prps %f',ps,pr,prps));  fflush(stdout);
 cmdst_neu = sprintf('-n %d -mat %s -pr %.16e -ps %.16e -scee %.16e', ...
                     p, matname, pr, ps, scee);
 id_stv = 0;
for savetimeinterval = s_stv
	id_stv = id_stv + 1;
	if id_stv == 1                 % output the first data, since it's longer
		cmdst_simu= sprintf('-t %.16e --save-interval %.16e -o "%s" --save-spike-interval %s --save-spike %s', ...
				     simu_time, savetimeinterval, output_name, output_ISI_name, output_RAS_name);
	else
		cmdst_simu= sprintf('-t %.16e --save-interval %.16e -o "%s"', ...
				     simu_time, savetimeinterval, output_name);
	end
	system([static_param, ' ',cmdst_neu, ' ',cmdst_simu]);

	% read volt data, binary format
%	fid = fopen(output_name, 'r');
%	X = fread(fid, [p, Inf], 'double');
%	fclose(fid);
%	[p, len] = size(X);

	ras = load('-ascii', output_RAS_name);
	for neuron_id = 1:p
		st = SpikeTrain(ras, len, neuron_id, 1, savetimeinterval);
		X(neuron_id,:) = st;
	end

	[oGC, oDe, R] = AnalyseSeries(X, s_od);

	prps_ps_stv_oGC(:,:,:, id_prps, id_ps, id_stv) = oGC;
	prps_ps_stv_oDe(:,:,:, id_prps, id_ps, id_stv) = oDe;
	prps_ps_stv_R  (:,:, id_prps, id_ps, id_stv) = R;
end  % stv
	prps_ps_aveISI(:, id_prps, id_ps) = load('-ascii', output_ISI_name);
	ISI_dis = zeros(p,length(hist_div));
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

%exit
