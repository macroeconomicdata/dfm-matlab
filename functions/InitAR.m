function [A_out, H, Q_out, R] = InitAR(Y, A, H, Q, R, frq, isdiff)

%InitCondAR()   Calculates initial estimates for autoregressive error
%               terms and shocks to the observation and transition equation
%
%Inputs:
%  Y: Standardized data
%  A: Transition matrix from model without AR errors
%  H: Loadings for model without AR errors
%  Q: Shocks to the transition equation for model without AR errors
%  R: Shocks to the observation equation for model without AR errors
%  frq: frequency mix
%  isdiff: logical, is the data differenced

%Output:
%  - A_out:   Transition matrix with AR parameters
%  - H:   Observation matrix (unchanged from input)
%  - Q:   Covariance for transition equation residuals
%  - R:   Covariance for observation equation residuals

[k,T] = size(Y);
sA = size(A,1);
p = sA/size(H,2);
HJ = get_HJ(H,frq,isdiff,p);
Zsmooth = runKF(Y, A, HJ, Q, R);

yBal = zeros(T,k);
for j=1:k
    yBal(:,j) = spline_fill_centered(Y(j,:));
end

%Shocks to observations
E = yBal - (HJ*Zsmooth(:,2:end))';
% AR components of transition matrix
a = zeros(1,k);
for j = 1:k
    a(j) = E(2:T,j)'*E(1:(T-1),j)/(E(1:T-1,j)'*E(1:(T-1),j));
end
A_out = zeros(sA+k,sA+k);
A_out(1:sA,1:sA) = A;
A_out(sA+1:sA+k, sA+1:sA+k) = diag(a); % AR error terms
E = E(2:T,:) - repmat(a,T-1,1).*E(1:T-1,:); % idiosyncratic errors
R = mean(E.^2); % variance of idio errors
Q_out = zeros(sA+k,sA+k);
Q_out(1:sA,1:sA) = Q;
Q_out(sA+1:sA+k, sA+1:sA+k) = diag(R);
return
end