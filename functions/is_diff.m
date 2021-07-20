function [isdiff]  = is_diff(Transformation)
  % Set the frequency mix
  % Input:
  %   String frequency: y, q, m, bw, w, d
  %
  % Output:
  %   Number of high frequency periods in each low frequency period to make
  %   the helper matrix J
  k = size(Transformation,1); %number of series
  isdiff = zeros(k,1);
  isdiff(strcmp('chg', Transformation)) = 1;
  isdiff(strcmp('pch', Transformation)) = 1;
  isdiff(strcmp('pca', Transformation)) = 1;
  isdiff = logical(isdiff);
  return
end
      