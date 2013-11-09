// mex isFileLocked.cpp

// input: filename to check
// output:
//  -1: file is not exist
//   0: file is exist and not locked
//   1: file is exist and locked

// see edit([matlabroot '/extern/examples/refbook/revord.c']);

#include <stdio.h>   // for fileno()
#include <unistd.h>  // for lockf()
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  if (nrhs != 1) {
    mexErrMsgTxt("no input file name.\n\
  usage: rt = isFileLocked(filename)\n\
  return -1 if file not exist;\n\
          0 if file exist and not locked;\n\
          1 file exist and locked.\n");
  }
  if (mxIsChar(prhs[0]) != 1 || mxGetM(prhs[0]) != 1) {
    mexErrMsgTxt("rt = isFileLocked(filename)\n  Input a string please.");
  }
  char *fn = mxArrayToString(prhs[0]);           // file name
  if (fn == NULL) {
    mexErrMsgTxt("isFileLocked: Can not convert input to string.");
  }
  FILE *f_tmp = fopen(fn,"r");
  double rt = 0;
  if (f_tmp == NULL) {
    rt = -1;                                     // file does not exist
  } else {
    if (lockf(fileno(f_tmp),F_TEST,0) == -1) {
      rt = 1;                                    // file is locked
    }
    fclose(f_tmp);
  }
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *y = mxGetPr(plhs[0]);
  *y = rt;

  return;
}
