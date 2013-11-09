
function flushstdout()
  persistent x;
  if isempty(x)
    x = exist('OCTAVE_VERSION','builtin');
  end
  if x
    fflush(stdout);
  end

