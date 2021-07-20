
%% PROCEDURES -------------------------------------------------------------
% Note: Kalman filter (runKF()) is in the 'functions' folder

function  [H_new, R_new, A_new, Q_new, loglik] = EMstepAR(Y, A, H, Q, R, p, frq, isdiff, blocks)
%EMstep    Applies EM algorithm for parameter reestimation
%
%  Syntax:
%    [C_new, R_new, A_new, Q_new, Z_0, V_0, loglik]
%    = EMstep(y, A, C, Q, R, Z_0, V_0, r, p, R_mat, q, nQ, i_idio, blocks)
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
%    Z_0:    Initial values of factors
%    V_0:    Initial value of factor covariance matrix
%    m:      Number of common factors for each block (e.g. vector [1 1 1 1])
%    p:      Number of lags in transition equation
%    blocks: Block structure for each series (i.e. for a series, the structure
%            [1 0 0 1] indicates loadings on the first and fourth factors)
%
%  Output:
%    H_new: Updated observation matrix
%    R_new: Updated covariance matrix for residuals of observation matrix
%    A_new: Updated transition matrix
%    Q_new: Updated covariance matrix for residuals for transition matrix
%    Z_0:   Initial value of state
%    V_0:   Initial value of covariance matrix
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
% Store series/model values
[k, T] = size(Y);
m = size(H,2);
sA = size(A,1); 
sA = sA - k; % less the AR error terms

%% ESTIMATION STEP: Compute the (expected) sufficient statistics for a single Kalman filter sequence

% Running the Kalman filter and smoother with current parameters
% Note that log-liklihood is NOT re-estimated after the runKF step: This
% effectively gives the previous iteration's log-likelihood
% For more information on output, see runKF
HJ = get_HJ(H, frq, isdiff, p); 
HH = [HJ, eye(k)];
[Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(Y, A, HH, Q, R);
% Vsmooth gives the variance of contemporaneous factors
% VVsmooth gives the covariance of factors at one lag for Watson Engle
% adjustments

%% Normalize

scl = mean((Zsmooth(1:m,:)).^2, 2); % scaling makes sure factors do not explode
scl = diag((scl).^-.5);
sscl = kron(eye(p),scl);
Scl = kron(eye(sA/m),scl);
Zsmooth(1:sA,:) = Scl*Zsmooth(1:sA,:);

% scl = 1/mean(mean(Zsmooth(1:m,:).^2)); % don't normalize AR error term factors
% Zsmooth(1:sA,:) = Zsmooth(1:sA,:)*sqrt(scl);

%% MAXIMIZATION STEP (TRANSITION EQUATION)
% See (Banbura & Modugno, 2010) for details.

%%% 2A. UPDATE FACTOR PARAMETERS  ----------------------------
% Initialize output
A_new = A;
Q_new = Q;
mp = m*p;
% ESTIMATE FACTOR PORTION OF Q, A
% Note: EZZ, EZZ_BB, EZZ_FB are parts of equations 6 and 8 in BM 2010

% E[f_t*f_t' | Omega_T]
EZZ = Zsmooth(1:mp, 2:end) * Zsmooth(1:mp, 2:end)'...
    + sscl*sum(Vsmooth(1:mp, 1:mp, 2:end) ,3)*sscl; % WE adjustment

% E[f_{t-1}*f_{t-1}' | Omega_T]
EZZ_BB = Zsmooth(1:mp, 1:end-1)*Zsmooth(1:mp, 1:end-1)'...
        + sscl*sum(Vsmooth(1:mp, 1:mp, 1:end-1), 3)*sscl; % WE adjustment

% E[f_t*f_{t-1}' | Omega_T]
EZZ_FB = Zsmooth(1:m, 2:end)*Zsmooth(1:mp, 1:end-1)'...
    + scl*sum(VVsmooth(1:m, 1:mp, :), 3)*sscl; % WE adjustment

% Equation 6: Estimate VAR(p) for factor
A_new(1:m,1:mp) = EZZ_FB/EZZ_BB; % VAR coeficients

% Equation 8: Covariance matrix of residuals of VAR
Q_new(1:m,1:m) = (EZZ(1:m,1:m) - A_new(1:m,1:mp)* EZZ_FB') / T; %shortcut --- very clever (again)

%% Update AR(1) components
% Below: estimate the idiosyncratic component (for eqns 6, 8 BM 2010)
% Just a bunch of AR(1) regressions

% sA + j locates where the AR component is for series j

% Everything from Vsmooth and VVsmooth are Watson Engle (1983) adjustment
% terms... otherwise this would just be OLS

for j = 1:k
    EZZ = sum(Zsmooth(sA+j,2:end).^2) + sum(Vsmooth(sA+j, sA+j, 2:end), 3);
    EZZ_BB = sum(Zsmooth(sA+j,1:end-1).^2) + sum(Vsmooth(sA+j, sA+j, 1:end-1), 3);
    EZZ_FB = Zsmooth(sA+j, 2:end)*Zsmooth(sA+j, 1:end-1)' + ...
        sum(VVsmooth(sA+j, sA+j, :), 3);
    A_new(sA+j, sA+j) = EZZ_FB/EZZ_BB;  % Equation 6
    Q_new(sA+j,sA+j) = (EZZ - A_new(sA+j,sA+j)* EZZ_FB) / T;
end

%% 3 MAXIMIZATION STEP (observation equation)

%%% INITIALIZATION AND SETUP ----------------------------------------------

% LOADINGS
H_new = H;

for j = 1:k % Loop through observables
    fq = frq(j);
    y = Y(j,:) - Zsmooth(sA+j,2:end);
    y_idx = ~isnan(y);
    y_obs = y(y_idx);
    lblock = logical(blocks(j,:));
    if fq==1
        Z_obs = Zsmooth(1:m,2:end); %drop pre sample value Z_0
        Z_obs = Z_obs(lblock,y_idx); %Z_obs where y observed
        V_obs = sum(Vsmooth(1:m,1:m,logical([0,y_idx])),3); %Vsmooth where y observed
        V_obs = scl(lblock, lblock)*V_obs(lblock, lblock)*scl(lblock, lblock); % for zero restrictions
    else
        J = helper_mat(fq,isdiff(j),m,sA);
        Z_obs = J*Zsmooth(1:sA,2:end);
        Z_obs = Z_obs(lblock,y_idx);
        V_obs = J*sum(Vsmooth(1:sA,1:sA,logical([0,y_idx])),3)*J';
        V_obs = scl(lblock, lblock)*V_obs(lblock, lblock)*scl(lblock, lblock);
    end
    V_ar = sum(Vsmooth(sA+j,sA+j,logical([0,y_idx])),3); %WE adjustment for AR(1) error term
    EZZ = Z_obs*Z_obs' + V_obs;
    h = (y_obs*Z_obs')/ EZZ;
    H_new(j,lblock) = h;
    R(j) = ((y_obs-h*Z_obs)*(y_obs-h*Z_obs)' + h*V_obs*h' + V_ar)/size(y_obs,2);
end
R_new = R;

return
end