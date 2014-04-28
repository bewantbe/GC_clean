% get structure from name of the network

function [network, matname] = getnetwork(netstr)

pathdir = fileparts(mfilename('fullpath'));
if isempty(netstr)
    matname = '-';
    network = [1];
else
    matname = [pathdir, '/network/', netstr, '.txt'];
    network = load('-ascii', matname);
end
