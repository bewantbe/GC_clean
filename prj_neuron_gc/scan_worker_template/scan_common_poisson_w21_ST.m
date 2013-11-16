% scan: percentage of common poisson input, fix others
% will save GC analyse results to .mat
tic();

signature = 'data_scan/comm_w21_ST';     % to distinguish different parallel program instances (also dir)

% simulation type
mode_IF   = 'IF';
mode_ST   = 1;
mode_comm = 1;

% simulation parameters
netstr = 'net_2_2';
scee   = 0.01;
pr     = 8.0;
ps     = 0.005;
stv    = 0.5;
dt     = 1.0/32;
iseed  = 1;                % set random seed
simu_time = 1e7;

% scan range
s_common_poisson = linspace(1e-15, 1-1e-15, 50);  % common poisson input

% analyse range
maxod  = 99;
s_od   = 1:maxod;          % order calculate range
hist_div = 0:0.5:200;      % ISI statistics range

rand('state', iseed);

if mode_comm == 1
    static_param = '../raster_tuning_co -ng -q --bin-save --RC-filter';
else
    static_param = '../raster_tuning -ng -q --bin-save --RC-filter';
end
data_path = ['data/', signature, '_'];

matname = ['network/', netstr, '.txt'];
neu_network = load('-ascii', matname);
p = size(neu_network, 1);
output_name     = [data_path, 'volt_', netstr, '.dat'];
output_ISI_name = [data_path, 'ISI_',  netstr, '.txt'];
output_RAS_name = [data_path, 'RAS_',  netstr, '.txt'];
cmdst_neu = sprintf('-n %d -mat %s -pr %.16e -ps %.16e -scee %.16e', ...
                    p, matname, pr, ps, scee);
cmdst_simu= sprintf('-t %.16e --save-interval %.16e -dt %.17e -o "%s" --save-spike-interval %s --save-spike %s', ...
                    simu_time, stv, dt, output_name, output_ISI_name, output_RAS_name);

prps = pr*ps;
len = round(simu_time/stv);
fprintf('ps %f, pr %f, prps %f\n',ps,pr,prps);  fflush(stdout);

s_aveISI  = zeros(p,                   length(s_common_poisson));
s_ISI_dis = zeros(p, length(hist_div), length(s_common_poisson));
s_oGC     = zeros(p, p, length(s_od),  length(s_common_poisson));
s_oDe     = zeros(p, p, length(s_od),  length(s_common_poisson));
s_R       = zeros(p, p*(maxod+1),      length(s_common_poisson));
for id_s_comm = 1:length(s_common_poisson)
    comm = s_common_poisson(id_s_comm);
disp(['comm: ', num2str(comm)]);  fflush(stdout);
    extst = sprintf('-seed %d --pr-mul %s', int32(2.0^31*rand(1)), num2str([(1-comm)*ones(1,p),comm],'%.16e '));
    system([static_param, ' ',cmdst_neu, ' ',cmdst_simu,' ',extst]);

    if mode_ST == 0
        fid = fopen(output_name, 'r');
        X = fread(fid, [p, Inf], 'double');
        fclose(fid);
        [p, len] = size(X);
    else
        ras = load('-ascii', output_RAS_name);
        for neuron_id = 1:p
            st = SpikeTrain(ras, len, neuron_id, 1, stv);
            X(neuron_id,:) = st;
        end
    end

    [oGC, oDe, R] = AnalyseSeries(X, s_od);

    aveISI = load('-ascii', output_ISI_name);
    ISI_dis = zeros(p,length(hist_div));
    ras = load('-ascii', output_RAS_name);
    for kk = 1:p
        tm = ras(ras(:,1)==kk, 2);
        tm = tm(2:end) - tm(1:end-1);
        ISI_dis(kk,:) = hist(tm, hist_div);
    end

    s_aveISI (:,    id_s_comm) = aveISI;
    s_ISI_dis(:,:,  id_s_comm) = ISI_dis;
    s_oGC    (:,:,:,id_s_comm) = oGC;
    s_oDe    (:,:,:,id_s_comm) = oDe;
    s_R      (:,:,  id_s_comm) = R;
end

datamatname = sprintf('%s_%s_sc=%g_pr=%.2f_ps=%g_t=%.1e_stv=%.1e.mat', signature, netstr, scee, pr, ps, simu_time, stv);
save('-v6', datamatname,...
     's_aveISI', 's_ISI_dis', 's_oGC', 's_oDe', 's_R', 's_common_poisson', 's_od', 'hist_div',...
     'maxod', 'dt', 'stv', 'len', 'simu_time', 'ps', 'pr', 'scee', 'netstr', 'p', 'iseed',...
     'mode_ST', 'mode_IF', 'mode_comm', 'signature');

toc();

cmdst = sprintf('notify-send "work \\\"%s\\\" done."', signature);
system(cmdst);

%exit
