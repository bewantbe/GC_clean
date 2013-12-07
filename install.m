% Compile all Matlab mex files used in this project
% Some files relie on Eigen (http://eigen.tuxfamily.org, we are using version 3.1.4). So you need to put the headers in correct place.

GC_CAL_HOME = fileparts(mfilename('fullpath'));

is_octave = exist('OCTAVE_VERSION','builtin');
is_linux  = isunix();  % ignore Mac

if is_octave
  flushstdout = @() fflush(stdout);
else
  flushstdout = @() 0;
end

url_download_eigen = 'http://bitbucket.org/eigen/eigen/get/3.1.4.zip';
dir_lib_ext = [GC_CAL_HOME,'/extern_lib/'];
if ~exist([dir_lib_ext,'eigen3'], 'file')
  fprintf('downloading and unziping Eigen...');  flushstdout();
  unzip(url_download_eigen,'./extern_lib/');
  movefile([dir_lib_base,'eigen-eigen-36bf2ceaf8f5'],[dir_lib_ext,'eigen3']);
  fprintf('done\n');  flushstdout();
end

disp('Compiling...');  flushstdout();
cd([GC_CAL_HOME,'/GCcal']);
if is_octave
  common_cmd_prefix = 'CXXFLAGS="-O3 -march=native -fopenmp -std=c++11"  LDFLAGS="-march=native -fopenmp" mkoctfile ';
  system([common_cmd_prefix,'--mex -o gendata_linear gendata_linear_v4_mex.cpp']);
  system([common_cmd_prefix,'--mex getcovzpdhded.cpp']);
  system([common_cmd_prefix,'MAfilter_v5.cpp']);
else
  common_cmd_prefix = ['mex CXXFLAGS="\$CXXFLAGS -std=c++11 -I',dir_lib_ext,'" CXXOPTIMFLAGS="-O3 -march=native" '];
  eval([common_cmd_prefix, 'gendata_linear_v4_mex.cpp -output gendata_linear']);
  eval([common_cmd_prefix, 'getcovzpdhded.cpp']);
  eval([common_cmd_prefix, 'MAfilter_v5.cpp']);
end

cd([GC_CAL_HOME,'/GCcal_spectrum']);
mex GaussianConvGrid.cpp
mex lowest_smooth_number_exact.cpp
mex lowest_smooth_number_fast.cpp

if is_octave
  cd([GC_CAL_HOME,'/prj_neuron_gc']);
  mkoctfile hist2dnn.cpp
  mex isFileLocked.cpp
end

% tip: use
%   find -path './.hg' -prune -o \( -name '*.oct' -o -name '*.mex*' \) -print
% to find all compiled files
