% scan: fix ps

% scan every thing (length - network - scee - ps - pr*ps - stv)
% save some simple analyse results to .mat
t0 = tic();
tocs = @(st) fprintf('%s: t = %6.3fs\n', st, toc());

ext_T = 1e4;  % avoid head effect

if isempty(strfind(upper(signature),upper('expIF')))
  dt_std = 1.0/32;
  model_neu = 'IF';
else
  dt_std = 0.004;
  model_neu = 'ExpIF';
end
mode_ST = isempty(strfind(upper(signature),upper('ST')))==0;  % using Spike Train?
if ~exist('extst','var')
    extst = '';
end
if strcmpi(model_neu,'ExpIF')
    extst = ['ExpIF ',extst];
end

signature0 = signature;
datamatname = signature;
if (length(s_net)==1)
  datamatname = sprintf('%s_%s', datamatname, s_net{1});
end
if (length(s_scee)==1)
  datamatname = sprintf('%s_sc=%g', datamatname, s_scee(1));
end
if (length(s_time)==1)
  datamatname = sprintf('%s_t=%.1e', datamatname, s_time(1));
end
%datamatname = [signature, '_info.mat'];  % the old way
datamatname = [datamatname, '_info.mat'];
save('-v7', datamatname, 'signature0', 's_net', 's_time', 's_scee', 's_prps', 's_ps', 's_stv', 's_od', 'hist_div', 'maxod', 'extst', 'dt_std');

data_dir_prefix = ['data/', signature, '_'];
ccnt = 0;

for net_id = 1:length(s_net)
 netstr = s_net{net_id};
 neu_network = getnetwork(netstr);
 p = size(neu_network, 1);
for simu_time = s_time
for id_scee=1:length(s_scee)
 scee = s_scee(id_scee);
 prps_ps_stv_oGC = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_oDe = zeros(p, p, length(s_od), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_stv_R   = zeros(p, p*(maxod+1), length(s_prps), length(s_ps), length(s_stv));
 prps_ps_aveISI  = zeros(p, length(s_prps), length(s_ps));
 prps_ps_ISI_dis = zeros(p, length(hist_div), length(s_prps), length(s_ps));

% prepare parallel
s_b_finished = false(length(s_stv), length(s_ps), length(s_prps));
s_b_launched = false(length(s_stv), length(s_ps), length(s_prps));
if ~exist('ncpu','var')
  [~, ncpu] = system('nproc');         % get number of cpus
  ncpu = max([str2num(ncpu), 1]);      % or leave one for other job?
end
while any(~s_b_finished(:))
 id_parallel = 0;
 goto_to_while_loop = false;

for id_prps=1:length(s_prps)
 prps = s_prps(id_prps);
for id_ps=1:length(s_ps)
 ps = s_ps(id_ps);
 pr = prps / ps;
for id_stv=1:length(s_stv)
  stv = s_stv(id_stv);
  len = round(simu_time/stv);

  id_parallel = id_parallel + 1;
  if s_b_finished(id_parallel)
    continue
  else
    if s_b_launched(id_parallel)
      [X, ISI, ras] = gendata_neu(netstr, scee, pr, ps, simu_time+ext_T, stv, ['read ' extst], data_dir_prefix);
      X = X(:, round(ext_T/stv)+1:end);  % cut head
      if ~isempty(ras)
        ras(ras(:,2)<=ext_T, :) = [];
        ras(:,2) = ras(:,2) - ext_T;
      end
      if isempty(X)
        pause(0.1);   % save CPU time
        continue
      end
    else
      % enough available cpu resources?
      if sum(s_b_launched(:))-sum(s_b_finished(:))>=ncpu 
        pause(0.1);   % save CPU time
        goto_to_while_loop = true;
        break;        % restart the loop, wait for the launched ones
      end
      gendata_neu(netstr, scee, pr, ps, simu_time+ext_T, stv, ['new ', extst, ' &'], data_dir_prefix);
      s_b_launched(id_parallel) = true;
      fprintf('id_parallel=%d launched.\n', id_parallel);  flushstdout();
      continue
    end
  end

  if mode_ST
    for neuron_id = 1:p
      st = SpikeTrain(ras, len, neuron_id, 1, stv);
      X(neuron_id,:) = st;
    end
  end

  [oGC, oDe, R] = AnalyseSeriesFast(X, s_od);

  prps_ps_stv_oGC(:,:,:, id_prps, id_ps, id_stv) = oGC;
  prps_ps_stv_oDe(:,:,:, id_prps, id_ps, id_stv) = oDe;
  prps_ps_stv_R  (:,:, id_prps, id_ps, id_stv) = R;
  if id_stv==1
    for id_p = 1:p  % recount ISI
      ISI(id_p) = simu_time/(sum(ras(:,1)==id_p,1));
    end
    prps_ps_aveISI(:, id_prps, id_ps) = ISI;
    ISI_dis = zeros(p,length(hist_div));
    for kk = 1:p
      tm = ras(ras(:,1)==kk, 2);
      tm = tm(2:end) - tm(1:end-1);
      ISI_dis(kk,:) = hist(tm, hist_div);
    end
    prps_ps_ISI_dis(:,:,id_prps, id_ps) = ISI_dis;
  end
  s_b_finished(id_parallel) = true;
  % save disk space
  gendata_neu(netstr, scee, pr, ps, simu_time+ext_T, stv, ['rm ' extst], data_dir_prefix);
  fprintf('id_parallel=%d finished.\n', id_parallel);  flushstdout();
  ccnt = ccnt + 1;
  fprintf('  -- ps %f, pr %f, prps %f, cnt=%d\n',ps,pr,prps,ccnt); flushstdout(); 
end  % for stv
  if goto_to_while_loop
    break;
  end
end  % for ps
  if goto_to_while_loop
    break;
  end
end  % for prps

end % while parallel

 datamatname = sprintf('%s_%s_sc=%g_t=%.1e.mat', signature, netstr, scee, simu_time);
 save('-v7', datamatname, 'prps_ps_stv_oGC', 'prps_ps_stv_oDe', 'prps_ps_stv_R', 'prps_ps_aveISI', 'prps_ps_ISI_dis');
end  % scee
end  % simu_time
end  % net

fprintf('Elapsed time is %6.3f\n', (double(tic()) - double(t0))*1e-6 );

%cmdst = sprintf('notify-send "work \\\"%s\\\" done."', signature);
%system(cmdst);

%exit
% vim: set ts=4 sw=4 ss=4
