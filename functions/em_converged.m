%--------------------------------------------------------------------------

function [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_decreased)
%em_converged    checks whether EM algorithm has converged
%
%  Syntax:
%    [converged, decrease] = em_converged(loglik, previous_loglik, threshold, check_increased)
%
%  Description:
%    em_converged() checks whether EM has converged. Convergence occurs if
%    the slope of the log-likelihood function falls below 'threshold'(i.e.
%    f(t) - f(t-1)| / avg < threshold) where avg = (|f(t)| + |f(t-1)|)/2
%    and f(t) is log lik at iteration t. 'threshold' defaults to 1e-4.
%
%    This stopping criterion is from Numerical Recipes in C (pg. 423).
%    With MAP estimation (using priors), the likelihood can decrease
%    even if the mode of the posterior increases.
%
%  Input arguments:
%    loglik: Log-likelihood from current EM iteration
%    previous_loglik: Log-likelihood from previous EM iteration
%    threshold: Convergence threshhold. The default is 1e-4.
%    check_decreased: Returns text output if log-likelihood decreases.
%
%  Output:
%    converged (numeric): Returns 1 if convergence criteria satisfied, and 0 otherwise.
%    decrease (numeric): Returns 1 if loglikelihood decreased.

%% Instantiate variables

% Threshhold arguments: Checks default behavior
if nargin < 3, threshold = 1e-4; end
if nargin < 4, check_decreased = 1; end

% Initialize output
converged = 0;
decrease = 0;

%% Check if log-likelihood decreases (optional)

if check_decreased
    if loglik - previous_loglik < -1e-3 % allow for a little imprecision
        fprintf(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik);
        decrease = 1;
    end
end

%% Check convergence criteria

delta_loglik = abs(loglik - previous_loglik);  % Difference in loglik
avg_loglik = (abs(loglik) + abs(previous_loglik) + eps)/2;

if (delta_loglik / avg_loglik) < threshold,
    converged = 1;  % Check convergence
end

end