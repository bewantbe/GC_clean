% Get the adjacency matrix from its (file) name
% Assume the matrix is stored in plain text file
% You may specify the dir path in `pathdir' if it's not in working dir and
%   not in dir of this file (./network)
%
% If the input `netstr' is adjacency matrix, this function will return its
% hash as file name in `matname'. You may specify its dir in pathdir also.
% network = netstr in this case.

function [network, matname] = getnetwork(netstr, pathdir)

e = filesep;
if ~exist('pathdir','var')
  pathdir = '.';    % search current dir
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

if isnumeric(netstr)
  % so netstr is already the adjacency matrix ?
  network = netstr;
  if size(network,1) ~= size(network,2) || ndims(network) ~= 2
    error('Not a square matrix!');
  end
  % give this matrix a name
  % the program can than save the matrix there
  matname = [pathdir, 'net_', num2str(length(network)),...
             '_0X', BKDRHash(mat2str(network))];
else
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

end
