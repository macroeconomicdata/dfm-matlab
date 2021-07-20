% Simulate data in DFM form

% Observation equation
% Y_t = H*F_t + e_t

% Transition equation
% F_t = A F_{t-1} + u_t

function [Y, X] = sim_dfm_mixed_frq(T)
B = [.4, -.2 , .2, .2, .1, 0; .3, .3 .1 .2, .1, .1];
A = comp_form(B);
% eig(A);
Z = zeros(T, 6); % 2 factors x 3 lags
q = [1, .5; .5, .8];
E = randn(T,2)*q;
H = [2, -1; 0, 2; .3, 1; -1, 1; 1, .5; .5, 1; 1.5, .1; 1, 3; .2, .5; 1, 1];
for t=2:T
    Z(t,:) = Z(t-1,:)*A' + [E(t,:), zeros(1,4)];
end
frq = [3;3;3;3;1;1;1;1;1;1]; % first 4 series are quarterly, otherwise monthly
is_diff = false(10,1); % no differenced data
HJ = get_HJ(H, frq, is_diff, 3); % p = 3 for 3 lags
Y = Z*HJ' + randn(T, 10);
idx = repmat([0;1;1], ceil(T/3),1); % index of NaN values for low frequency data
Y(logical(idx(1:T)),1:4) = NaN; % drop intermediate low frequency values
X = Z(:,1:2);
end
