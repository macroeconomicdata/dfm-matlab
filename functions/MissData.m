function [y,C,R,L]  = MissData(y,C,R)
% Syntax:
% Description:
%   Eliminates the rows in y & matrices C, R that correspond to missing 
%   data (NaN) in y
%
% Input:
%   y: Vector of observations at time t
%   C: Observation matrix
%   R: Covariance for observation matrix residuals
%
% Output:
%   y: Vector of observations at time t (reduced)     
%   C: Observation matrix (reduced)     
%   R: Covariance for observation matrix residuals
%   L: Used to restore standard dimensions(n x #) where # is the nr of 
%      available data in y
  
  % Returns 1 for nonmissing series
  ix = ~isnan(y);
  
  % Index for columns with nonmissing variables
  e  = eye(size(y,1));
  L  = e(:,ix);

  % Removes missing series
  y  = y(ix);
  
  % Removes missing series from observation matrix
  C  =  C(ix,:);  
  
  % Removes missing series from transition matrix
  R  =  R(ix,ix);

end

