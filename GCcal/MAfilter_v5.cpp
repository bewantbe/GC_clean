// compile with
// CXXFLAGS='-O3 -march=native -fopenmp -std=c++11'  LDFLAGS='-march=native -fopenmp' mkoctfile MAfilter_v5.cpp
#include <octave/oct.h>

#define EIGEN_NO_DEBUG           // then there is no range checking in Eigen
#define CMAKE_BUILD_TYPE Release
#include <eigen3/Eigen/Dense>

#define HELP_STR "MA filter (B, X)"
DEFUN_DLD (MAfilter_v5, args, nargout, HELP_STR)
{
  int nargin = args.length();
  if (nargin!=2) {
    print_usage();
    return octave_value_list();
  }

  const Matrix B = args(0).matrix_value();
  const Matrix X = args(1).matrix_value();
  const int n = X.dim2();
  const int p = B.dim1();
  const int pm= B.dim2();
  const int m = (int)floor((double)pm/p+0.5);
  if (X.dim1() != B.dim1()) {
    error("MAfilter: B is p * (p*m) matrix, X is p * n matrix.");
    return octave_value_list();
  }

  // copy MA coefficients, in transposed inverse-time order
  Eigen::MatrixXd B_roll(pm, p);
  for (int j=0; j<m; j++) {
    for (int l=0; l<p*p; l++) {
      B_roll((m-1-j)*p+l/p, l%p) = B(l%p, j*p+l/p);
    }
  }

  const Eigen::Map<const Eigen::RowVectorXd, Eigen::Aligned> X_map(X.data(), p*n);
  Matrix e_X(p, n);  // final result
  Eigen::Map<Eigen::RowVectorXd, Eigen::Aligned> Y_map(e_X.fortran_vec(), p*n);

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

  return octave_value(e_X);
}
