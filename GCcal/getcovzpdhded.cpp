// copy from prj_GC_clean/test_code
// mex getcovpdhded.cpp
//In octave
// CXXFLAGS='-O3 -march=native -fopenmp -std=c++11'  LDFLAGS='-march=native -fopenmp' mkoctfile --mex getcovpdhded.cpp
#include <mex.h>

#define EIGEN_NO_DEBUG
#include <eigen3/Eigen/Dense>
// use Col Major

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 4 || ! mxIsNumeric(prhs[3])) {
    mexErrMsgTxt("getcovpdhded(X, t_bg, t_ed, m)");
  }
  size_t p    = (size_t)mxGetM(prhs[0]);  // mrows
  size_t len  = (size_t)mxGetN(prhs[0]);  // ncols
  size_t t_bg = (size_t)mxGetScalar(prhs[1]);
  size_t t_ed = (size_t)mxGetScalar(prhs[2]);
  int m       =    (int)mxGetScalar(prhs[3]);
  if (m+1 > len) {
    mexErrMsgTxt("m too large!");
  }
  if (t_bg < m+1) {
    mexErrMsgTxt("t_bg too small!");
  }
  if (t_ed > len-m) {
    mexErrMsgTxt("t_ed too large!");
  }

  size_t m1   = m + 1;
  int len_ed = len-t_ed;
  int t_bg_m = t_bg-1-m;
  plhs[0] = mxCreateDoubleMatrix(m1*p, m1*p, mxREAL);
  Eigen::Map<Eigen::MatrixXd> covz(mxGetPr(plhs[0]), m1*p, m1*p);
  const Eigen::Map<const Eigen::MatrixXd> X(mxGetPr(prhs[0]),p,len);
  for (int i1=0; i1<=m; i1++) {
    for (int i2=i1; i2<=m; i2++) {
      covz.block(i1*p,i2*p,p,p) =
        X.block(0,m-i1, p,t_bg_m+i1) * 
        X.block(0,m-i2, p,t_bg_m+i1).transpose() +
        X.block(0,t_ed      , p,len_ed-i1) *
        X.block(0,t_ed-i2+i1, p,len_ed-i1).transpose();
      if (i1!=i2)
        covz.block(i2*p,i1*p,p,p) = covz.block(i1*p,i2*p,p,p).transpose();
    }
  }
}
