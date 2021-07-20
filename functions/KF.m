%______________________________________________________________________
function S = KF(Y, A, HJ, Q, R)
%  Applies fast Kalman filter (see Dornbush Koopman (2012))
%
%  Syntax:
%    S = SKF(Y, A, C, Q, R)
%
%  Description:
%    Kfilter() applies the Kalman filter for the approximate model, meaning
%    that R is a vector, not a matrix.

%  Input parameters:
%    Y: k-by-nobs matrix of input data
%    A: m-by-m transition matrix 
%    C: k-by-m observation matrix
%    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
%    R: k-by-1 vector of variances for shocks to observations(e_t)
%
%  Output parameters:
%    S.Zm: m-by-nobs matrix, prior/predicted factor state vector
%          (S.Zm(:,t) = Z_t|t-1)
%    S.ZmU: m-by-(nobs+1) matrix, posterior/updated state vector
%           (S.Zm(t+1) = Z_t|t)
%    S.Vm: m-by-m-by-nobs array, prior/predicted covariance of factor
%          state vector (S.Vm(:,:,t) = V_t|t-1)  
%    S.VmU: m-by-m-by-(nobs+1) array, posterior/updated covariance of
%           factor state vector (S.VmU(:,:,t+1) = V_t|t)
%    S.loglik: scalar, value of likelihood function
%    S.k_t: m-by-m output for smoothing initial (nobs) obs
  
%% INITIALIZE OUTPUT VALUES ---------------------------------------------
  % Output structure & dimensions of state space matrix
  sA = size(A,2); %size of transition matrix (it's square)
  
  % Outputs size for data matrix. "number of observations"
  [k,nobs]  = size(Y);
  
  % Instantiate output
  S.Zm  = nan(sA, nobs);       % Z_t | t-1 (prior)
  S.Vm  = nan(sA, sA, nobs);    % V_t | t-1 (prior)
  S.ZmU = nan(sA, nobs+1);     % Z_t | t (posterior/updated)
  S.VmU = nan(sA, sA, nobs+1);  % V_t | t (posterior/updated)
  S.loglik = 0; % initialize log likelihood
  S.UD = zeros(sA,k,nobs); % where we store factor updates

%% SET INITIAL VALUES ----------------------------------------------------
  % Initial values follow Hamilton (1994)
  Z = zeros(sA, 1);  % Z_0|0 (In below loop, Zu gives Z_t | t)
  V = long_run_var(A, Q);  % V_0|0 (In below loop, Vu guvse V_t | t)
  
  % Store initial values
  S.ZmU(:,1)    = Z; % store presample factors
  S.VmU(:,:,1)  = V; % store presample variance

%% KALMAN FILTER PROCEDURE ----------------------------------------------
  for t = 1:nobs
      %%% CALCULATING PRIOR DISTIBUTION----------------------------------
      
      % Use transition eqn to create prior estimate for factor
      % i.e. Z = Z_t|t-1
      Z   = A * Z; % prediction for factors in next period
      
      % Prior covariance matrix of Z (i.e. V = V_t|t-1)
      %   Var(Z) = Var(A*Z + u_t) = Var(A*Z) + Var(\epsilon) = 
      %   A*Vu*A' + Q
      V   = A * V* A' + Q; % variance of factor prediction
      V   =  0.5 * (V+V');  % Trick to make symmetric
      
      % Store covariance and observation values for t-1 (priors)
      S.Zm(:,t)   = Z; % store predicted factors
      S.Vm(:,:,t) = V; % store variance of prediction
      
      %%% CALCULATING POSTERIOR DISTRIBUTION ----------------------------
       
      y_t = Y(:,t); % extract observatoins in the current period
      ix = find(~isnan(y_t)); % observed values in y_t
      % Check if y_t contains no data. If so, replace Zu and Vu with prior.
      if ~isempty(y_t)
          for jj=1:size(ix,1)
              j = ix(jj); % data which is in fact observed (index of)
              y_j = y_t(j); % scalar data
              h_j = HJ(j,:); % corresponding row of the H (loadings) matrix
              s = h_j*V*h_j' + R(j); % varaince of forecast for y_j
              pe = y_j - h_j*Z; % forecast error, forecast is h_j*Z
              K = V*h_j'/s; % no matrix inversions!! Kalman gain for this specific series
              Z = Z + K*pe; % update estimated factors based on observatoin j
              V = V - K*h_j*V; % update variance of factors based on observation j
              S.loglik = S.loglik - (log(2*pi) + log(s) + pe*pe/s)/2; % likelihood of the single observation
              S.UD(:,j,t) = K * pe; %store update contributions for each series
          end % once we loop through all observatoins we break out and use the updated values of Z and V
      end
      
      %%% STORE OUTPUT----------------------------------------------------

      % Store covariance and state values for t (posteriors)
      % i.e. Zu = Z_t|t   & Vu = V_t|t
      S.ZmU(:,t+1)    = Z; % store these variables for all observatoins in period t
      S.VmU(:,:,t+1)  = V;
  end 
  % Store Kalman gain for last period k_t for smoothing
  ix = ~isnan(Y(:,nobs));
  if all(~ix)
      S.k_t = zeros(sA,sA);
  else
      H_t = HJ(ix,:);
      C  = V * H_t';  
      P   = H_t * C + diag(R(ix));
      K = (P\C')';
      S.k_t = K * H_t;
  end
  
end