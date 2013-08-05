% Give spike train from spike time data

function X = SpikeTrains(ras, p, len, stv)
if ~exist('stv', 'var') || isempty(stv)
  stv = 0.5;
end

for neuron_id = 1:p
  st = SpikeTrain(ras, len, neuron_id, 1, stv);
  X(neuron_id,:) = st;
end
