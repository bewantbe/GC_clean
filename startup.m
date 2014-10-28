% Startup file for user-defined options

GC_CAL_HOME = fileparts(mfilename('fullpath'));

addpath(GC_CAL_HOME);
addpath([GC_CAL_HOME,'/GCcal']);
addpath([GC_CAL_HOME,'/GCcal_spectrum']);
addpath([GC_CAL_HOME,'/tools_and_utilities']);
addpath([GC_CAL_HOME,'/prj_neuron_gc']);
addpath([GC_CAL_HOME,'/prj_neuron_gc/scan_worker_template']);

disp('GC calculation package path added.');
disp(' ');

% check whether the mex files are up to date
cppfiles={
'GCcal/gendata_linear_v4_mex.cpp'
'GCcal/getcovzpdhded.cpp'
'GCcal/MAfilter_v5_mex.cpp'
'GCcal_spectrum/GaussianConvGrid.cpp'
'GCcal_spectrum/lowest_smooth_number_exact.cpp'
'GCcal_spectrum/lowest_smooth_number_fast.cpp'
};
mexfiles={
'GCcal/gendata_linear'
'GCcal/getcovzpdhded'
'GCcal/MAfilter'
'GCcal_spectrum/GaussianConvGrid'
'GCcal_spectrum/lowest_smooth_number_exact'
'GCcal_spectrum/lowest_smooth_number_fast'
};
for id_f=1:length(cppfiles)
  path_cpp = [GC_CAL_HOME '/' cppfiles{id_f}];
  path_mex = [GC_CAL_HOME '/' mexfiles{id_f} '.' mexext];
  dcpp=dir(path_cpp);
  dmex=dir(path_mex);
  if isempty(dcpp)
    disp(['GC package: "' path_cpp '" not found!']);
    error('GC package: your copy of this package is broken.');
    break;
  end
  if isempty(dmex)
    if exist('OCTAVE_VERSION','builtin')
      path_mex = [GC_CAL_HOME '/' mexfiles{id_f} '.oct'];
      dmex=dir(path_mex);
    end
    if isempty(dmex)
      disp(['GC package: "' path_mex '" not found!']);
      warning('GC package: run `install'' to generate all necessary mex files.');
    end
    break;
  end
  if dcpp.datenum > dmex.datenum
    warning('GC package: some mex file(s) are not up to date, run `install'' to re-generate them.');
    break;
  end
end

clear('cppfiles', 'dmex', 'mexfiles', 'path_mex', 'path_cpp', 'id_f', 'dcpp', 'ans');
