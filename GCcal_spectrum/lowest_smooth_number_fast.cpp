// find the approximate minimum number larger than input
// that only have factors 2,3,5
// borrow from code/numerical/small_factor_number

#define COMPILE_MODE 2

typedef unsigned int UINT;
const UINT MAX_UINT = 4294967295;

const UINT tbmul[] = {32, 36, 40, 45, 48, 50, 54, 60};
const UINT smtab[] = {2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32};
const UINT smtab_size = 18;  // =sizeof(smtab)/sizeof(UINT)

UINT SmallFactorApp5(UINT n)
{
  UINT n2=1;
  while (n2<n) {
    n2 <<= 1;
  }
  if (n2==n)
    return n;

  UINT sm = n2;
  if (n > smtab[smtab_size-1]) {
    n2 >>= 6;
    for (int id=0; id<smtab_size; id++) {
      UINT m = n2 * tbmul[id];
      if (m >= n) {
        sm = m;
        break;
      }
    }
  } else {
    int id = -1;
    while (smtab[++id]<n) {}
    sm = smtab[id];
  }

  return sm;
}

#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 1 || ! mxIsNumeric(prhs[0])) {
    mexErrMsgTxt("expects a number");
  }

  double v = *mxGetPr(prhs[0]);
  if (!(v>0 && v<MAX_UINT/2)) {
    mexErrMsgTxt("expects a number within range [1, 2147483647]");
  }
  UINT n = (UINT)(v+0.5);

  mxArray *ov = mxCreateDoubleMatrix (1, 1, mxREAL);
  double *data = mxGetPr (ov);
  *data = SmallFactorApp5(n);
  plhs[0] = ov;
}
