% parallel generate test

mode_IF = 'IF';
mode_ST = 0;
netstr = 'net_2_2';
scee = 0.01;
pr = 2;
ps = 0.006;
simu_time = 1e5;
extst = '--RC-filter -q';

if ~exist('extst','var')
    extst = '';
end
if strcmpi(mode_IF,'ExpIF')
    extst = ['ExpIF ',extst];
end
if exist('dt','var')
    extst = [extst, sprintf(' -dt %.16e',dt)];
end
if ~exist('stv','var')
    stv = 1/2;
end

od = 20;
tic;

network = getnetwork(netstr);
p = size(network,1);
prps = pr*ps;
s_ps = linspace(0.001, 0.02, 100);
s_p_gc = zeros(p*p-p, length(s_ps));
s_b_finished = false(size(s_ps));
s_b_launched = false(size(s_ps));
[~, ncpu] = system('nproc');         % get number of cpus
ncpu = max([str2num(ncpu), 1]);      % or leave one for other job?
% delete the data files in this loop
for id_ps=1:length(s_ps)
  ps = s_ps(id_ps);
  pr = prps/ps;
  gendata_neu(netstr, scee, pr, ps, simu_time, stv, ['rm ' extst]);
end
% parallel scan!
while any(~s_b_finished(:))
for id_ps=1:length(s_ps)
  ps = s_ps(id_ps);
  pr = prps/ps;
  %fprintf('   ps=%e, ps=%e\n', pr, ps);
  if s_b_finished(id_ps)
    continue
  else
    if s_b_launched(id_ps)
      %fprintf('id_ps=%d try read data.\n', id_ps);  flushstdout();
      [X, ISI, ras] = gendata_neu(netstr, scee, pr, ps, simu_time, stv, ['read ' extst]);
      if isempty(X)
        %fprintf('id_ps=%d fail to read data.\n', id_ps);  flushstdout();
        pause(0.1);   % save CPU time
        continue
      end
    else
      % enough available cpu resources?
      if sum(s_b_launched)-sum(s_b_finished)>=ncpu 
        %disp('lack of resources');
        pause(0.1);   % save CPU time
        break;        % restart the loop, wait for the launched ones
      end
      gendata_neu(netstr, scee, pr, ps, simu_time, stv, [extst, ' &']);
      s_b_launched(id_ps) = true;
      fprintf('id_ps=%d launched.\n', id_ps);  flushstdout();
      continue
    end
  end
  gc = nGrangerT(X, od);
  gc(eye(size(gc))==1) = [];
  s_p_gc(:, id_ps) = gc;
  s_b_finished(id_ps) = true;
  % save disk space
  gendata_neu(netstr, scee, pr, ps, simu_time, stv, ['rm ' extst]);
  fprintf('id_ps=%d finished.\n', id_ps);  flushstdout();
end %for
end %while any

toc
