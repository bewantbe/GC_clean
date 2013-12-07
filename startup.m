% Startup file for user-defined options

GC_CAL_HOME = fileparts(mfilename('fullpath'));

addpath(GC_CAL_HOME);
addpath([GC_CAL_HOME,'/GCcal']);
addpath([GC_CAL_HOME,'/GCcal_spectrum']);
addpath([GC_CAL_HOME,'/tools_and_utilities']);
addpath([GC_CAL_HOME,'/prj_neuron_gc']);
addpath([GC_CAL_HOME,'/prj_neuron_gc/scan_worker_template']);

disp('GC calculation packages loaded.');
disp(' ');
