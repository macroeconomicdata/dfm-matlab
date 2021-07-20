% Simulate data in DFM form

% Observation equation
% Y_t = H*F_t + e_t

% Transition equation
% F_t = A F_{t-1} + u_t

function [Y, X] = sim_dfm(T)
B = [.4, -.2 , .2, .2, .1, 0; .3, .3 .1 .2, .1, .1];
% A = comp_form(B);
% eig(A);
X = zeros(T, 2);
q = [1, .5; .5, .8];
E = randn(T,2)*chol(q);
H = [2, -1; 0, 2; .3, 1; -1, 1; 1, .5; .5, 1; 1.5, .1; 1, 3; .2, .5; 1, 1];
for t=4:T
    X(t,:) = [X(t-1,:), X(t-2,:), X(t-3,:)]*B' + E(t,:);
end
Y = X*H' + randn(T, 10);
end
