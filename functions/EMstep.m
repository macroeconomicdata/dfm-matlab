
function  [H_new, R_new, A_new, Q_new, loglik] = EMstep(Y, A, H, Q, R, p, frq, isdiff, blocks)
%EMstep    Applies EM algorithm for parameter reestimation
%
%
%  Description:
%    EMstep reestimates parameters based on the Estimation Maximization (EM)
%    algorithm. This is a two-step procedure:
%    (1) E-step: the expectation of the log-likelihood is calculated using
%        previous parameter estimates.
%    (2) M-step: Parameters are re-estimated through the maximisation of
%        the log-likelihood (maximize result from (1)).
%
%    See "Maximum likelihood estimation of factor models on data sets with
%    arbitrary pattern of missing data" for details about parameter
%    derivation (Banbura & Modugno, 2010). This procedure is in much the
%    same spirit.
%
%  Input:
%    y:      Series data
%    A:      Transition matrix
%    H:      Observation matrix
%    Q:      Covariance for transition equation residuals
%    R:      Covariance for observation matrix residuals
%    p:      Number of lags in transition equation
%    frq:    frequency of each series in y
%    isdiff: logical (T/F) indicating if series in y is differenced
%    
%
%  Output:
%    H_new: Updated observation matrix
%    R_new: Updated covariance matrix for residuals of observation matrix
%    A_new: Updated transition matrix
%    Q_new: Updated covariance matrix for residuals for transition matrix
%    loglik: Log likelihood
%
% References:
%   "Maximum likelihood estimation of factor models on data sets with
%   arbitrary pattern of missing data" by Banbura & Modugno (2010).
%   Abbreviated as BM2010
%
%

%% Initialize preliminary values

% Store series/model values
[k, T] = size(Y); % size of the data
[~,m] = size(H); % number of factors
sA = size(A,1); % total number of stacked factors (size of A)
pp = sA/m; % total number of lags in the model

%% ESTIMATION STEP: Compute the (expected) sufficient statistics for a single Kalman filter sequence

% Running the Kalman filter and smoother with current parameters
% Note that log-liklihood is NOT re-estimated after the runKF step: This
% effectively gives the previous iteration's log-likelihood
% For more information on output, see runKF
HJ =get_HJ(H, frq, isdiff, p); 

% A = Spec.A;
% HJ = Spec.HJ;
% Q = Spec.Q;
% R = Spec.R;

[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(Y, A, HJ, Q, R); %diag(R)
% Vsmooth gives the variance of contemporaneous factors
% VVsmooth gives the covariance of factors at one lag for Watson Engle
% adjustments

%% Normalize

scl = mean((Zsmooth(1:m,:)).^2, 2); % scaling makes sure factors do not explode
scl = diag((scl).^-.5);
sscl = kron(eye(p),scl);
Scl = kron(eye(pp),scl);
Zsmooth = Scl*Zsmooth;


%% MAXIMIZATION STEP (TRANSITION EQUATION)
% See (Banbura & Modugno, 2010) for details.

%%% 2A. UPDATE FACTOR PARAMETERS  ----------------------------
% Initialize output
A_new = A; % only replace B
Q_new = Q; % only replace q (i.e. covariance for contemporaneous shocks)
mp = m*p; % size of B
% ESTIMATE FACTOR PORTION OF Q, A
% Note: EZZ, EZZ_BB, EZZ_FB are parts of equations 6 and 8 in BM 2010

% E[Z_t*Z_t' | Omega_T]
% Variance of contemporaneous factors P in the notes
EZZ = Zsmooth(1:mp, 2:end) * Zsmooth(1:mp, 2:end)'...
    + sscl*sum(Vsmooth(1:mp, 1:mp, 2:end) ,3)*sscl; % WE adjustment
% scl is in there because we didn't scale things prior to running the
% filter/smoother

EZZ = (EZZ + EZZ')/2; % just to avoid rounding error

% E[Z_{t-1}*Z_{t-1}' | Omega_T]
% P_{t-1|T} in the notes
EZZ_BB = Zsmooth(1:mp, 1:end-1)*Zsmooth(1:mp, 1:end-1)'...
        + sscl*sum(Vsmooth(1:mp, 1:mp, 1:end-1), 3)*sscl; % WE adjustment
    
EZZ_BB = (EZZ_BB + EZZ_BB')/2;

% E[X_t*Z_{t-1}' | Omega_T]
% sum of C_{t|T} in the notes
EZZ_FB = Zsmooth(1:m, 2:end)*Zsmooth(1:mp, 1:end-1)'...
    + scl*sum(VVsmooth(1:m, 1:mp, :), 3)*sscl; % WE adjustment

% Equation 6: Estimate VAR(p) for factor
A_new(1:m,1:mp) = EZZ_FB/EZZ_BB; % VAR coeficients

% Equation 8: Covariance matrix of residuals of VAR
% same as slide 50 in dfm slides
q = (EZZ(1:m,1:m) - A_new(1:m,1:mp)* EZZ_FB') / T; %shortcut --- very clever (again)
Q_new(1:m,1:m) =  (q + q')/2;

%% 3 MAXIMIZATION STEP (observation equation)

%%% INITIALIZATION AND SETUP ----------------------------------------------

% LOADINGS
H_new = H;

for j = 1:k % Loop through observables
    fq = frq(j); % get the frequency of this series
    y = Y(j,:); % the actual observed data
    y_idx = ~isnan(y); % find values that are observed
    y_obs = y(y_idx);
    lblock = logical(blocks(j,:)); % which factors can this series load on?
    if fq==1
        Z_obs = Zsmooth(1:m,2:end); %drop pre sample value Z_0, columns are time in y
        Z_obs = Z_obs(lblock,y_idx); %Z_obs where y observed for those factors y can load on
        % sum of P^x_{t|T} in the slides
        V_obs = sum(Vsmooth(1:m,1:m,logical([0,y_idx])),3); %Vsmooth where y observed
        V_obs = scl(lblock,lblock)*V_obs(lblock, lblock)*scl(lblock,lblock); % for zero restrictions
    else
        J = helper_mat(fq,isdiff(j),m,sA);
        Z_obs = J*Zsmooth(1:sA,2:end);
        Z_obs = Z_obs(lblock,y_idx);
        V_obs = J*sum(Vsmooth(1:sA,1:sA,logical([0,y_idx])),3)*J';
        V_obs = scl(lblock,lblock)*V_obs(lblock, lblock)*scl(lblock,lblock); % for zero restrictions
    end
    EZZ = Z_obs*Z_obs' + V_obs;
    EZZ = (EZZ + EZZ')/2; % get rid of rounding error
    h = (y_obs*Z_obs')/ EZZ; % almost OLS, just adjusting for the fact that factors are estimated not observed
    H_new(j,lblock) = h; % plug estimated values of h in
    R(j) = (sum((y_obs-h*Z_obs).^2)+ h*V_obs*h')/size(y_obs,2); % shocks to this series in observation equation
end
R_new = R;

return
end