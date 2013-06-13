% Startup file for user-defined options

ph = fileparts(mfilename('fullpath'));

addpath([ph,'/GCcal']);
addpath([ph,'/experimental_tools']);
addpath([ph,'/prj_neuron_gc']);
addpath([ph,'/SpectrumCal']);
addpath([ph,'/npGCcal']);

disp('GC calculation packages loaded.');
disp(' ');
