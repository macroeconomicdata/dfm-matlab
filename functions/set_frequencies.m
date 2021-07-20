function [frq]  = set_frequencies(string_frq)
  % Set the frequency mix
  % Input:
  %   String frequency: y, q, m, bw, w, d
  %
  % Output:
  %   Number of high frequency periods in each low frequency period to make
  %   the helper matrix J
  if iscell(string_frq)
      string_frq = cellfun(@(x) x(1), string_frq); % the function dfm only wants the first letter
  end
  k = size(string_frq,1); %number of series
  frq = 28*ones(k,1); %if nothing is found, default is monthly
  % numerical value for each frequency
  frq(ismember(string_frq, 'd')) = 1;
  frq(ismember(string_frq, 'w')) = 7;
  frq(ismember(string_frq, 'b')) = 14; % biweekly
  frq(ismember(string_frq, 'm')) = 28;
  frq(ismember(string_frq, 'q')) = 84;
  frq(ismember(string_frq, 'y')) = 336;
  min_frq = min(frq);
  frq = frq/min_frq; % number of high frequency preiods in a low frequency period
  %finalizing
  if min_frq == 1
      frq(frq == 28) = 31;
      frq(frq == 84) = 91;
      frq(frq == 336) = 365;
  elseif min_frq == 7
      frq(frq == 12) = 13;
      frq(frq == 48) = 52;
  elseif min_frq == 14
      frq(frq == 24) = 26;
  end
  return
end
      