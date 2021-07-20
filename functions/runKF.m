
function [Zsmooth, Vsmooth, VVsmooth, loglik, Update] = runKF(Y, A, HJ, Q, R)
%runKF()    Applies Kalman filter and fixed-interval smoother
%
%  Syntax:
%    [zsmooth, Vsmooth, VVsmooth, loglik] = runKF(Y, A, C, Q, R, Z_0, V_0)
%
%  Description:
%    runKF() applies a Kalman filter and fixed-interval smoother. The
%    script uses the following model:
%           Y_t = C_t Z_t + e_t for e_t ~ N(0, R)
%           Z_t = A Z_{t-1} + mu_t for mu_t ~ N(0, Q)

%  Throughout this file:
%    'm' denotes the number of elements in the state vector Z_t.
%    'k' denotes the number of elements (observed variables) in Y_t.
%    'nobs' denotes the number of time periods for which data are observed.
%
%  Input parameters:
%    Y: k-by-nobs matrix of input data
%    A: m-by-m transition matrix 
%    C: k-by-m observation matrix
%    Q: m-by-m covariance matrix for transition equation residuals (mu_t)
%    R: k-by-k covariance for observation matrix residuals (e_t)
%    Z_0: 1-by-m vector, initial value of state
%    V_0: m-by-m matrix, initial value of state covariance matrix
%
%  Output parameters:
%    zsmooth: k-by-(nobs+1) matrix, smoothed factor estimates
%             (i.e. zsmooth(:,t+1) = Z_t|T)
%    Vsmooth: k-by-k-by-(nobs+1) array, smoothed factor covariance matrices
%             (i.e. Vsmooth(:,:,t+1) = Cov(Z_t|T))
%    VVsmooth: k-by-k-by-nobs array, lag 1 factor covariance matrices
%              (i.e. Cov(Z_t,Z_t-1|T))
%    loglik: scalar, log-likelihood
%
%  References:
%  - QuantEcon's "A First Look at the Kalman Filter"
%  - Adapted from replication files for:
%    "Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia 
%    Reichlin), in Michael P. Clements and David F. Hendry, editors, Oxford 
%    Handbook on Economic Forecasting.
%
% The software can be freely used in applications. 
% Users are kindly requested to add acknowledgements to published work and 
% to cite the above reference in any resulting publications

S = KF(Y, A, HJ, Q, R);  % Kalman filter
S = Ksmooth(A, S);       % Fixed interval smoother

% Organize output 
Zsmooth = S.ZmT;
Vsmooth = S.VmT;
VVsmooth = S.VmT_1;
loglik = S.loglik;
Update = S.UD;

end