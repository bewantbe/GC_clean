% Get the adjacency matrix from its (file) name
% Assume the matrix is stored in plain text file
% You may specify the dir path in `pathdir' if it's not in working dir and
%   not in dir of this file (./network)

function [network, matname] = getnetwork(netstr, pathdir)

if ~ischar(netstr)
  error('Input should be the name of the matrix');
end

e = filesep;
if ~exist('pathdir','var')
  pathdir = '';    % search default dir (usually working dir)
end
% always consider pathdir as dir name
if ~isempty(pathdir) && pathdir(end) ~= '/' && pathdir(end) ~= e
  pathdir = [pathdir e];
end
if isempty(netstr)
  matname = '-';
  network = [1];
  return
end

pathdir0 = fileparts(mfilename('fullpath'));
% path candidates, priority from high to low
s_matname = {...
  [pathdir, netstr, '.txt'],...
  [pathdir, netstr],...
  [pathdir0, e, 'network', e, netstr, '.txt'],...
  [pathdir0, e, 'network', e, netstr]...
};
for matname = s_matname
  matname = matname{1};
  fprintf('f=%s\n', matname);
  if isempty(dir(matname))
    continue
  end
  network = load('-ascii', matname);
  break
end
if ~exist('network','var')
  error('Network file not found!');
end

end
