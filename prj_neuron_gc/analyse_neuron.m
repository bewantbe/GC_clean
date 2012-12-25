% Example Matlab M file
% Show the use of raster_tuning.exe and nGrangerT
% Run this script under nGranger/prj_neuron_gc/
% See readme(raster_tuning).txt for more information about 
% the program raster_tuning or use raster_tuning -h for a short help
% (run system('../raster_tuning -h'); in Matlab)

% add path to nGrangerT.m etc. so we can run nGrangerT here
addpath ../GCcal
tic();
% number of neurons to calculate(Only Ex. type)
num_neu_ex = 3;
% time length of simulation. there are 2 data for each millisecond
% In order to get a acceptable result, at least 3e4, 1e5 is better,
% and 1e6 will get an impress result
simu_time = 1e6;

% preparing a matrix for cortical network
matname = 'network/net_3_06.txt';
%mfcn = rand(num_neu_ex, num_neu_ex) > 0.5;
%save('-ascii', matname, 'mfcn');

scee = 0.02;    % cortical strength
pr = 0.8;       % poisson input rate
ps = 0.010;     % poisson input strength
% the command line to calculate neuron volt data,
cmdst = sprintf('./raster_tuning -ng -q --bin-save -n %d -t %g -mat %s -pr %.16e -ps %.16e -scee %.16e --save-spike-interval data/spike_int.txt', ...
num_neu_ex, simu_time, matname, pr, ps, scee);

% generate the data
system(cmdst);

% load the data (in text form)
%load('data/staffsave.txt');
%X = staffsave';

fid = fopen('data/staffsave.txt', 'r');
X = fread(fid, [num_neu_ex, Inf], 'double');
fclose(fid);
%Y=X;
%X = [mean(Y(1:3,:)); mean(Y(4:6,:))];

od = 10;
% calculate Granger Causality in time domain with order od
G1 = pos_nGrangerT2(X, od);            % unbiased(-cov) version
disp(G1);

% save the results
%save('-ascii', 'gc_result_1.txt', 'G1');


% analyse the results

% load true cortial matrix
ms = load('-ascii', matname);
a = G1(ms > 0.5 & eye(num_neu_ex) < 0.5);       % the GC of "one"
b = G1(ms < 0.5 & eye(num_neu_ex) < 0.5);       % the GC of "zero"
amax = max(a);
amin = min(a);
bmax = max(b);
bmin = min(b);
k1 = amin/bmax;
k2 = amax/amin;
k3 = bmax/bmin;
disp('compare the min "one" and max "zero"');
fprintf('%g / %g = %f\n', amin, bmax, k1); % compare the two results, demand >1, bigger is better

disp('compare the max "one" and min "one"');
fprintf('%g / %g = %f\n', amax, amin, k2);

disp('compare the max "zero" and min "zero"');
fprintf('%g / %g = %f\n', bmax, bmin, k3);

toc()
