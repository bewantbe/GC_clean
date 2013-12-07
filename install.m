% Compile all Matlab mex files used in this project
% You should have a compiler that support C++11

GC_CAL_HOME = fileparts(mfilename('fullpath'));

is_octave = exist('OCTAVE_VERSION','builtin');
is_linux  = isunix();  % ignore Mac

if is_octave
  flushstdout = @() fflush(stdout);
else
  flushstdout = @() 0;
end

url_download_eigen = 'https://bitbucket.org/eigen/eigen/get/3.1.4.zip';
dir_lib_ext = [GC_CAL_HOME,'/extern_lib/'];
if ~exist([dir_lib_ext,'eigen3'], 'file')
  fprintf('downloading and unziping Eigen...');  flushstdout();
  %unzip(url_download_eigen,'extern_lib');  % Seems Matlab does not follow the redirect hint (HTTP 301) in the url.
  urlwrite(url_download_eigen,'eigen-3.1.4.zip');
  unzip('eigen-3.1.4.zip','extern_lib');
  movefile([dir_lib_ext,'eigen-eigen-36bf2ceaf8f5'],[dir_lib_ext,'eigen3']);
  fprintf('done\n');  flushstdout();
end

disp('Compiling...');  flushstdout();
cd([GC_CAL_HOME,'/GCcal']);
if is_octave
  common_cmd_prefix = 'CXXFLAGS="-O3 -march=native -fopenmp -std=c++11"  LDFLAGS="-march=native -fopenmp" mkoctfile ';
  system([common_cmd_prefix,'--mex -o gendata_linear gendata_linear_v4_mex.cpp']);
  system([common_cmd_prefix,'--mex getcovzpdhded.cpp']);
  system([common_cmd_prefix,'-o MAfilter MAfilter_v5.cpp']);
else
  common_cmd_prefix = ['mex CXXFLAGS="\$CXXFLAGS -std=c++11 -I',dir_lib_ext,'" CXXOPTIMFLAGS="-O3 -march=native" '];
  eval([common_cmd_prefix, 'gendata_linear_v4_mex.cpp -output gendata_linear']);
  eval([common_cmd_prefix, 'getcovzpdhded.cpp']);
  eval([common_cmd_prefix, 'MAfilter_v5_mex.cpp -output MAfilter']);
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
