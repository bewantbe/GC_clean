// find the minimum number larger than input that only have factors 2,3,5
// borrow from code/numerical/small_factor_number

typedef unsigned int UINT;
const UINT MAX_UINT = 4294967295;
UINT smallest;

void SmallFactorApp2(UINT n, UINT m, UINT p)
{
  if (n>m) {
    switch (p) {
      case 2 : SmallFactorApp2(n, m*2, 2);
      case 3 : SmallFactorApp2(n, m*3, 3);
      case 5 : SmallFactorApp2(n, m*5, 5);
    }
  } else {
    if (m<smallest) smallest = m;
  }
}

#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  if (nrhs != 1 || ! mxIsNumeric(prhs[0])) {
    mexErrMsgTxt("expects a number");
  }

  double v = *mxGetPr(prhs[0]);
  if (!(v>0 && v<MAX_UINT/5)) {
    mexErrMsgTxt("expects a number within range [1, 858993459]");
  }
  UINT n = (UINT)(v+0.5);

  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  double *data = mxGetPr(plhs[0]);
  smallest = MAX_UINT;
  SmallFactorApp2(n, 1, 2);
  *data = smallest;
}

