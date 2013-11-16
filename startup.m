% Startup file for user-defined options

ph = fileparts(mfilename('fullpath'));

addpath(ph);
addpath([ph,'/GCcal']);
addpath([ph,'/GCcal_spectrum']);
addpath([ph,'/experimental_tools']);
addpath([ph,'/prj_neuron_gc']);
addpath([ph,'/prj_neuron_gc/scan_worker_template']);

disp('GC calculation packages loaded.');
disp(' ');
