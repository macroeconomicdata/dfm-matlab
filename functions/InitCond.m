function [A, H, Q, R] = InitCond(X,m,p,frq,isdiff, blocks)

%InitCond()      Calculates initial conditions for parameter estimation
%
%  Description:
%    Given standardized data and model information, InitCond() creates
%    initial parameter estimates. These are intial inputs in the EM
%    algorithm, which re-estimates these parameters using Kalman filtering
%    techniques.
%
%Inputs:
% X: Standardized data
% m: Number of common factors
% p: Number of lags in transition equation
% frq: frequency mix
% isdiff: logical, is the series differenced
% blocks: zero restrictions on loadings
%Output:
%  - A:   Transition matrix
%  - H:   Observation matrix
%  - Q:   Covariance for transition equation residuals
%  - R:   Covariance for observation equation residuals

% This is just estimation by principal components or "two-step" estimation

[T,k] = size(X);
xBal = zeros(T,k);
% fill in missing values via cubic spline
for j=1:k
    xBal(:,j) = spline_fill_centered(X(:,j));
end

H = zeros(k, m); % initialize loadings
xTemp = xBal; 
z = zeros(T,m); %initialize factors
for j = 1:m
    lblock = logical(blocks(:,j)); % loading restrictions for this factor
    [G, ~] = eigs(cov(xTemp(:,lblock)), 1, 'lm'); % Initial guess for loadings 
    H(lblock, j) = G; % G is the eigenvector which will be our loading 
    z(:,j) = xTemp(:,lblock)*G; % z is our factor 
    xTemp = xTemp - z(:,j)*H(:,j)'; % removing the variance explained by this factor from the data
end
% VAR for transition equation
Z = stack_obs(z,p,true);
sV = size(Z,2);
B = (z(p+1:T,:)'*Z)/(Z'*Z + eye(sV)); %transition matrix
E = z(p+1:T,:)-Z*B'; % shocks to factors
%Adjusting for differenced low frequency data
% -- This tells us how many lags we need in the case of mixed frequencies--
lags = frq;
lags(isdiff,:) = arrayfun(@(x)(2*x-1),frq(isdiff,:));
pp = max([lags;p]);
% ------------------------------------
sA = m*pp; %size of A matrix
Q = zeros(sA,sA); % dimensions must line up: shocks to transition equation
Q(1:m,1:m) = E'*E/(T-p);
%Shocks to observations
E = xBal - z*H';
E = E(2:T,:);
R = mean(E.^2) + 1; % shocks to observation equatoin 

if pp>p
    B = [B, zeros(m,m*(pp-p))]; % again make sure dimensions line up
end
A = zeros(sA);
A(1:m*pp, 1:m*pp) = comp_form(B); % DFM is specified for the companion form of the transition matrix

return
end