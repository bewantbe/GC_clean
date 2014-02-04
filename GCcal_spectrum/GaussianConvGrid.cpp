#include <mex.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// return non-negtive remainder
int mod(int a, int m)
{
  int b = a % m;
  return (b<0) ? m+b : b;
}

// Calculate Gaussian spread at time T[j+k*N] with strength X[j+k*N]
// j=0..N, k=0..p-1, there are p trials, length of each trial is N.
 
// Limitation:
//   Range of T should within [0, 2*pi]
//   subscripts range of X(N*p), T(N*p) and return value(M_r*p)
//     should smaller than MAX_INT (usually 2^31-1 = 2147483647)
//   M_sp  < sqrt(MAX_INT)
//   M_r   < sqrt(MAX_INT)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 5 || ! mxIsNumeric(prhs[0])) {
    mexErrMsgTxt("GaussianConvGrid(X, T, M_r, M_sp, tau)");
  }

  double *X = mxGetPr(prhs[0]);
  double *T = mxGetPr(prhs[1]);
  int M_r  = (int)mxGetScalar(prhs[2]);
  int M_sp = (int)mxGetScalar(prhs[3]);
  double tau = (double)mxGetScalar(prhs[4]);
  int N = mxGetM(prhs[0]);
  int p = mxGetN(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(M_r, p, mxREAL);
  double *f_r = mxGetPr(plhs[0]);

  double *exp3 = (double*)malloc((M_sp+1)*sizeof(double));
  if (exp3==NULL) {
    mexErrMsgTxt("GaussianConvGrid: no enough memory! (fail to malloc)");
  }
  for (int l=0; l<=M_sp; l++) {
    exp3[l] = exp(-l*l*(M_PI*M_PI)/(M_r*M_r)/tau);
  }
  for (int k=0; k<p; k++) {                   // variable loop
    for (int j=0; j<N; j++) {                 // source point loop
      int id = floor(T[j+k*N]/(2*M_PI)*M_r);
      double tj = T[j+k*N] - 2*M_PI*id/M_r;   // centralize the time
      double exp1 = exp(-tj*tj/4/tau);
      double exp20= exp(tj*M_PI/M_r/tau);
      double exp2 = 1;
      f_r[mod(id,M_r) + k*M_r] += exp1 * X[j+k*N];
      for (int l=1; l<=M_sp; l++) {           // spread Gaussian loop
        double tmp = exp1 * exp3[l] * X[j+k*N];
        exp2 *= exp20;
        f_r[mod(id+l,M_r) + k*M_r] += tmp * exp2;
        f_r[mod(id-l,M_r) + k*M_r] += tmp / exp2;
      }
    }
  }
  free(exp3);
}
