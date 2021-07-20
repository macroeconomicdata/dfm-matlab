%______________________________________________________________________
function S = Kfilter(Y, A, HJ, Q, R)
%  Applies Kalman filter
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
%    R: k-by-k matrix of variances for shocks to observations(e_t)
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
%    S.k_t: k-by-m Kalman gain
  
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
  S.loglik = 0;
  S.UD = zeros(sA,k,nobs);

%% SET INITIAL VALUES ----------------------------------------------------
  % Initial values follow Hamilton (1994)
  Zu = zeros(sA, 1);  % Z_0|0 (In below loop, Zu gives Z_t | t)
  Vu = long_run_var(A, Q);  % V_0|0 (In below loop, Vu guvse V_t | t)
  
  % Store initial values
  S.ZmU(:,1)    = Zu;
  S.VmU(:,:,1)  = Vu;

%% KALMAN FILTER PROCEDURE ----------------------------------------------
  for t = 1:nobs
      %%% CALCULATING PRIOR DISTIBUTION----------------------------------
      
      % Use transition eqn to create prior estimate for factor
      % i.e. Z = Z_t|t-1
      Z   = A * Zu;
      
      % Prior covariance matrix of Z (i.e. V = V_t|t-1)
      %   Var(Z) = Var(A*Z + u_t) = Var(A*Z) + Var(\epsilon) = 
      %   A*Vu*A' + Q
      V   = A * Vu* A' + Q; 
      V   =  0.5 * (V+V');  % Trick to make symmetric
      
      %%% CALCULATING POSTERIOR DISTRIBUTION ----------------------------
       
      % Removes missing series: These are removed from Y, C, and R
      yt = Y(:,t);
      ix = ~isnan(yt); % identify missing observations
      % Removes missing series
      yt  = yt(ix);  
      % Removes missing series from observation matrix
      H_t  =  HJ(ix,:);    
      % Removes missing series from transition matrix
      R_t  =  R(ix,ix);

      % Check if y_t contains no data. If so, replace Zu and Vu with prior.
      if isempty(yt)
          Zu = Z;
          Vu = V;
      else  
          % Steps for variance and population regression coefficients:
          % Var(c_t*Z_t + e_t) = c_t Var(A) c_t' + Var(u) = c_t*V *c_t' + R
          C  = V * H_t';  
          P   = H_t * C + R_t; %variance of observables
          Pi  = pinv(P); 
          
          % Matrix of population regression coefficients (QuantEcon eqn #4)
          K = C*Pi;  %Kalman gain K

          % Gives difference between actual and predicted observation
          % matrix values
          innov  = yt - H_t*Z;
          
          % Update estimate of factor values (posterior)
          Zu  = Z  + K * innov;
          
          % Update covariance matrix (posterior) for time t
          Vu  = V  - K * C';
          Vu   =  0.5 * (Vu+Vu'); % Approximation trick to make symmetric
          
          % Update log likelihood 
          S.loglik = S.loglik - (size(yt,1)*log(2*pi) + log(det(P)) + innov'*Pi*innov)/2;
          
          % Store Update Contribution for contemporanious factors
          S.UD(:,ix,t) = K .* kron(ones(sA,1), innov');
      end
      
      %%% STORE OUTPUT----------------------------------------------------
      
      % Store covariance and observation values for t-1 (priors)
      S.Zm(:,t)   = Z;
      S.Vm(:,:,t) = V;

      % Store covariance and state values for t (posteriors)
      % i.e. Zu = Z_t|t   & Vu = V_t|t
      S.ZmU(:,t+1)    = Zu;
      S.VmU(:,:,t+1)  = Vu;
  end 
  % Store Kalman gain for last period k_t for smoothing
  if isempty(yt)
      S.k_t = zeros(sA,sA);
  else
      S.k_t = K * H_t;
  end
  
end