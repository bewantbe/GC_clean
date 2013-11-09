#include <octave/oct.h>

#define HELP_STRING "Usage: nn = hist2dnn(x, y, bounds, nbins)"

DEFUN_DLD (hist2dnn, args, nargout, HELP_STRING)
{
  int nargin = args.length ();
  if (nargin != 4) {
    octave_stdout << HELP_STRING;
    return octave_value_list ();
  }
  double x_min = args(2).matrix_value()(0);
  double x_max = args(2).matrix_value()(1);
  int nbins = args(3).row_vector_value()(0);

  double k = (nbins-1)/(x_max-x_min);
  double b = x_min;
  double v0, v1;
  int id0, id1;

  Matrix n2d(nbins, nbins);
  for (int i=0; i<nbins; i++)
    for (int j=0; j<nbins; j++)
      n2d(i,j) = 0; 

  if (!args(0).is_real_matrix()) {
    octave_stdout << "non-compatible input type\n";
    return octave_value_list ();
  }
  Matrix X = args(0).matrix_value();
  Matrix Y = args(1).matrix_value();
  int len = X.dim1()>X.dim2() ? X.dim1() : X.dim2();
  for (octave_idx_type i=0; i<len; i++) {
    v0 = X(i);
    v1 = Y(i);

    if (v0>x_max)
      id0 = nbins-1;
    else if (v0<x_min)
      id0 = 0;
    else
      id0 = round(k*(v0 - b));

    if (v1>x_max)
      id1 = nbins-1;
    else if (v1<x_min)
      id1 = 0;
    else
      id1 = round(k*(v1 - b));

    n2d(id0, id1)++;
  }

  octave_value_list retval;
  retval(0) = n2d;
  return retval;
}
