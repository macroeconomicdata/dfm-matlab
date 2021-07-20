
function [Zsmooth, Vsmooth, VVsmooth, loglik, Update] = runKF(Y, A, HJ, Q, R)
%runKF()    Applies Kalman filter and fixed-interval smoother
%%  Description:
%    runKF() applies a Kalman filter and fixed-interval smoother. The
%    script uses the following model:
%           Y_t = C_t Z_t + e_t for e_t ~ N(0, R)
%           Z_t = A Z_{t-1} + mu_t for mu_t ~ N(0, Q)
%
%%  Input parameters:
%    Y: k-by-nobs matrix of input data
%    A: m-by-m transition matrix 
%    HJ: k-by-m observation matrix
%    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
%    R: k-by-k covariance for observation matrix residuals (e_t)
%
%  Output parameters:
%    Zsmooth: k-by-(nobs+1) matrix, smoothed factor estimates
%             (i.e. zsmooth(:,t+1) = Z_t|T)
%    Vsmooth: k-by-k-by-(nobs+1) array, smoothed factor covariance matrices
%             (i.e. Vsmooth(:,:,t+1) = Cov(Z_t|T))
%    VVsmooth: k-by-k-by-nobs array, lag 1 factor covariance matrices
%              (i.e. Cov(Z_t,Z_t-1|T))
%    loglik: scalar, log-likelihood
%    Update: Contribution of each observation to the update of each factor
%            for a period t

S = KF(Y, A, HJ, Q, R);  % Kalman filter
S = Ksmooth(A, S);       % Fixed interval smoother

% Organize output 
Zsmooth = S.ZmT;
Vsmooth = S.VmT;
VVsmooth = S.VmT_1;
loglik = S.loglik;
Update = S.UD;

end