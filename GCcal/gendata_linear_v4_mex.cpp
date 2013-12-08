// port gendata_linear_v4.cpp to matlab
// You should have Eigen installed
// Compile under Matlab:
//   Windows (Assume using MSVC as default compiler):
//     mex COMPFLAGS="$COMPFLAGS /fopenmp" LINKFLAGS="$LINKFLAGS /fopenmp" gendata_linear_v4_mex.cpp -output gendata_linear
//   UNIX-like:
//     mex CXXFLAGS="\$CXXFLAGS -std=c++11" CXXOPTIMFLAGS="-O3 -march=native -fopenmp"  LDFLAGS="\$LDFLAGS -fopenmp" gendata_linear_v4_mex.cpp -output gendata_linear
// Compile under Octave (at shell command window):
//   CXXFLAGS='-O3 -march=native -fopenmp -std=c++11' LDFLAGS='-march=native -fopenmp' mkoctfile --mex gendata_linear_v4_mex.cpp
#include <mex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// for random number
#include <random>
unsigned long GetRndSeed()
{
  std::random_device rd;
  return rd();
}

std::mt19937 eng;
std::normal_distribution<double> ran_gaussian;

inline double randn()
{
  return ran_gaussian(eng);
}

#define EIGEN_NO_DEBUG           // then there is no range checking in Eigen
#define CMAKE_BUILD_TYPE Release
#include <eigen3/Eigen/Dense>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs!=3 && nrhs!=4) {
    mexErrMsgTxt("usage: gendata_linear (A, De, n [, seed])\n");
  }

  const double *A = mxGetPr(prhs[0]);
  const double *D = mxGetPr(prhs[1]);
  const int n = (int)mxGetScalar(prhs[2]);
  const int p = mxGetM(prhs[0]);  // rows
  const int pm= mxGetN(prhs[0]);  // cols
  const int m = pm / p;

  plhs[0] = mxCreateDoubleMatrix(p, n, mxREAL);
  Eigen::Map< Eigen::RowVectorXd, Eigen::Aligned > X_map(mxGetPr(plhs[0]), p*n);

  Eigen::MatrixXd A_roll(pm, p);
  for (int j=0; j<m; j++) {
    for (int l=0; l<p*p; l++) {
      A_roll((m-1-j)*p+l/p, l%p) = A[l + p*p*j];
    }
  }

  unsigned long rng_seed = GetRndSeed();
  if (nrhs == 4) {
    rng_seed = (int)mxGetScalar(prhs[3]);
  }
  eng.seed(rng_seed);
  Eigen::RowVectorXd X_pre(2*pm);
  for (int j=0; j<pm; j++) {
    X_pre(j) = randn();
  }

  Eigen::MatrixXd hD(p,p);    // chol decomposition
  for (int j=0; j<p*p; j++) {
    hD(j) = D[j];
  }
  Eigen::LLT<Eigen::MatrixXd> lltOfD(hD);
  hD = lltOfD.matrixU();

  Eigen::RowVectorXd xt(p);
  int it=0;
  for (; it<pm; it+=p) {
    for (int j=0; j<p; j++) {
      xt(j) = randn();
    }
    X_pre.segment(it+pm, p).noalias() = xt*hD;
    X_pre.segment(it+pm, p).noalias()-= X_pre.segment(it, pm) * A_roll;
    X_map.segment(it, p) = X_pre.segment(it+pm, p);
  }

  const int np = n*p;
  for (; it<np; it+=p) {
    for (int j=0; j<p; j++) {
      xt(j) = randn();
    }
    X_map.segment(it, p).noalias() = xt*hD;
    X_map.segment(it, p).noalias()-= X_map.segment(it-pm, pm) * A_roll;
  }

}
