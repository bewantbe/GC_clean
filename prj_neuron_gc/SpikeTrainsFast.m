% Output spike train from spike time data

function X = SpikeTrainsFast(ras, p, len, stv, st_mode)
  if any(ras(:,1)>p)
    error('!! any(ras(:,1)>p)');
  end
  if ras(end, 2) >= stv*len || ras(1, 2) < 0
    % Assume ras(:,2) is sorted.
    error('!! ras(end, 2) >= stv*len || ras(1, 2) < 0');
  end
  if ~exist('st_mode', 'var') || st_mode == 0
    X = zeros(p, len);
    X( ras(:, 1) + p*(floor(ras(:, 2)/stv)) ) = 1;
  else
    % also count repeated values
    idv = ras(:, 1) + p*(floor(ras(:, 2)/stv));
    unv = unique(idv);
    X = zeros(p, len);
    X(unv) = histc(idv,unv);
  end
end

