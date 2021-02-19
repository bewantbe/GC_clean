% Compile all Matlab mex/Octave oct files used in this project
% You should have a compiler that support C++11
% Here assume you are using GCC or a GCC compatible compiler under unix-like system or using MSVC under MS Windows.

GC_CAL_HOME = fileparts(mfilename('fullpath'));

is_octave = exist('OCTAVE_VERSION','builtin');
is_linux  = isunix();  % ignore Mac

if is_octave
  flushstdout = @() fflush(stdout);
else
  flushstdout = @() 0;
end

url_download_eigen = 'https://gitlab.com/libeigen/eigen/-/archive/3.1.4/eigen-3.1.4.zip';
dir_lib_ext = [GC_CAL_HOME,filesep,'extern_lib',filesep];
if ~exist([dir_lib_ext,'eigen3'], 'file')
  fprintf('Downloading Eigen...');  flushstdout();
  urlwrite(url_download_eigen,'eigen-3.1.4.zip');  % Seems Matlab does not follow the redirect hint (HTTP 301) in the url.
  fprintf('Unziping Eigen...');  flushstdout();
  unzip('eigen-3.1.4.zip','extern_lib');
  % rename it, since I use #include <eigen3/***> in the cpp source code
  movefile([dir_lib_ext,'eigen-3.1.4'],[dir_lib_ext,'eigen3']);
  fprintf('done.\n');  flushstdout();
end

disp('Compiling...');  flushstdout();

if is_linux
    cd([GC_CAL_HOME,'/GCcal']);
    if is_octave
      common_cmd_prefix = 'CXXFLAGS="-O3 -march=native -fopenmp -std=c++11"  LDFLAGS="-march=native -fopenmp" mkoctfile ';
      system([common_cmd_prefix,'--mex -o gendata_linear gendata_linear_v4_mex.cpp']);
      system([common_cmd_prefix,'--mex getcovzpdhded.cpp']);
      system([common_cmd_prefix,'MAfilter_v5.cpp -o MAfilter']);
    else
      % Matlab
      common_cmd_prefix = ['mex CXXFLAGS="\$CXXFLAGS -std=c++11 -I',dir_lib_ext,'" CXXOPTIMFLAGS="-O3 -march=native -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" '];
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
else
    % In Windows. Assume you are using MSVC (>=2010)
    cd([GC_CAL_HOME,'/GCcal']);
    common_cmd_prefix = ['mex COMPFLAGS="$COMPFLAGS /I ',dir_lib_ext,'" OPTIMFLAGS="$OPTIMFLAGS /fopenmp" LINKFLAGS="$LINKFLAGS /fopenmp" '];
    eval([common_cmd_prefix, 'gendata_linear_v4_mex.cpp -output gendata_linear']);
    eval([common_cmd_prefix, 'getcovzpdhded.cpp']);
    eval([common_cmd_prefix, 'MAfilter_v5_mex.cpp -output MAfilter']);

    cd([GC_CAL_HOME,'/GCcal_spectrum']);
    mex GaussianConvGrid.cpp
    mex lowest_smooth_number_exact.cpp
    mex lowest_smooth_number_fast.cpp
end

cd(GC_CAL_HOME);

% tip(in unix-like): use
%   find -path './.hg' -prune -o \( -name '*.oct' -o -name '*.mex*' \) -print
% to find all compiled files
