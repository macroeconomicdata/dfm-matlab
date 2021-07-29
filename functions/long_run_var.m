function V_0 = long_run_var(A,Q)
% Calculate the long run variance of the transition equation given the
% paremeters A and Q
sA = size(A,1);
xx = eye(sA^2) - kron(A,A);
vQ = reshape(Q, (sA)^2, 1); % vectorize Q
V_0 = xx\vQ;
V_0 = reshape(V_0,sA,sA); 
V_0 = (V_0 + V_0')/2; %reduce rounding error
end

