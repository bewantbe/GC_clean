// port MAfilter_v5.cpp to matlab
// Compile under Matlab (you should have Eigen (eigen.tuxfamily.org) installed first):
//   Windows:
//   UNIX-like:
//     mex CXXFLAGS="\$CXXFLAGS -std=c++11" CXXOPTIMFLAGS="-O3 -march=native" MAfilter_v5_mex.cpp -output MAfilter
// Compile under Octave (at shell command window):
//   CXXFLAGS='-O3 -march=native -fopenmp -std=c++11'  LDFLAGS='-march=native -fopenmp' mkoctfile --mex MAfilter_v5.cpp -o MAfilter
#include <mex.h>

#define EIGEN_NO_DEBUG           // then there is no range checking in Eigen
#define CMAKE_BUILD_TYPE Release
#include <eigen3/Eigen/Dense>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs!=2) {
    mexErrMsgTxt("Multi-variable MA filter.\n  usage: MAfilter(B, X)\nSee also: filter\n");
  }

  const double *B = mxGetPr(prhs[0]);
  const double *X = mxGetPr(prhs[1]);
  const int n = mxGetN(prhs[1]);  // cols
  const int p = mxGetM(prhs[0]);  // rows
  const int pm= mxGetN(prhs[0]);  // cols
  const int m = (int)floor((double)pm/p+0.5);
  if (p != mxGetM(prhs[1])) {
    mexErrMsgTxt("MAfilter: B is p * (p*m) matrix, X is p * n matrix.\n");
  }

  // copy MA coefficients, in transposed inverse-time order
  Eigen::MatrixXd B_roll(pm, p);
  for (int j=0; j<m; j++) {
    for (int l=0; l<p*p; l++) {
      B_roll((m-1-j)*p+l/p, l%p) = B[l+p*p*j];
    }
  }

  const Eigen::Map<const Eigen::RowVectorXd, Eigen::Aligned> X_map(X, p*n);
  plhs[0] = mxCreateDoubleMatrix(p, n, mxREAL);
  Eigen::Map<Eigen::RowVectorXd, Eigen::Aligned> Y_map(mxGetPr(plhs[0]), p*n);

  // for head part
  Eigen::RowVectorXd X_pre(2*pm);
  X_pre.segment(0, pm).setZero();                // initial values
  X_pre.segment(pm, pm) = X_map.segment(0, pm);

  const int minmnp = (m<n ? m : n)*p;
  int it=0;
  for (; it<minmnp; it+=p) {
    Y_map.segment(it, p).noalias() = X_pre.segment(it+pm, p) + X_pre.segment(it, pm) * B_roll;
  }

  // remain part
  const int np = n*p;
  for (; it<np; it+=p) {
    Y_map.segment(it, p) = X_map.segment(it, p);
    Y_map.segment(it, p).noalias() += X_map.segment(it-pm, pm) * B_roll;
  }
}
