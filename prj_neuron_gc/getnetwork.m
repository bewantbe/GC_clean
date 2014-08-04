% Get the adjacency matrix from its (file) name.
% Assume the matrix is stored in plain text file. Which essentially can be
%   loaded  by `network = load('-ascii', [netstr, '.txt'])'
% It will search current working dir (or pathdir if specified) and
%   dir of this function.
% It also possible to specifies the full path in `netstr'.

function [network, matname] = getnetwork(netstr, pathdir)

if ~ischar(netstr)
  error('Input should be the name of the matrix');
end

if ~exist('pathdir','var')
  % Use default dir (usually working dir)
  pathdir = '';
end
% always consider pathdir as dir name
e = filesep;
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
  %fprintf('f=%s\n', matname);
  if ~exist(matname, 'file')
    continue
  end
  network = load('-ascii', matname);
  break
end
if ~exist('network','var')
  error('Network file not found!');
end

end
