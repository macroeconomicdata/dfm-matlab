function Res = dfm(X,X_pred,m,p,frq,isdiff,varargin) % blocks, threshold, ar_errors
%DFM()    Runs the dynamic factor model
%
%  Input arguments:
%    X: Transformed input data (i.e. log differenced where needed).
%    m: number of factors
%    p: number of lags
%    isdiff: logical (T/F) indicating if series in corresponding column is
%    differenced (for mixed frequency aggregation)
%
% Output Arguments:
%
%   Res - structure of model results with the following fields
%       .X_sm  Kalman-smoothed data where missing values are replaced by their expectation
%       .Z  Smoothed states. 
%       .H  Observation matrix. For low frequency data, C gives loadings for aggregated factors
%       .HJ Observation matrix times appropriate helper matrix J for each
%       series. This is what is actually used to get fitted values of
%       observables. 
%       .R: Covariance for observation matrix residuals
%       .A: Transition matrix
%       .Q: Covariance for transition equation residuals.
%       .Mx: Series mean
%       .Wx: Series standard deviation
%       .Z_0: Initial value of state
%       .V_0: Initial value of covariance matrix
%       .r: Number of common factors for each block
%       .p: Number of lags in transition equation
%
% References:
%
%   Marta Banbura, Domenico Giannone and Lucrezia Reichlin
%   Nowcasting (2010)
%   Michael P. Clements and David F. Hendry, editors,
%   Oxford Handbook on Economic Forecasting.

%% Store model parameters ------------------------------------------------

fprintf('Estimating the dynamic factor model (DFM) ... \n\n');

[T,k] = size(X); % get the size of the data
isdiff = logical(isdiff);

if nargin == 6
    blocks = ones(k,m);
    threshold = 1e-5;  % EM loop threshold (default value)
    ar_errors = false;
    varnames = (1:k)';
elseif nargin == 7
    blocks = varargin{1};
    threshold = 1e-5;  % EM loop threshold (default value)
    ar_errors = false;
    varnames = (1:k)';
elseif nargin == 8
    blocks = varargin{1};
    ar_errors = varargin{2};
    threshold = 1e-5;
    varnames = (1:k)';
elseif nargin == 9
    blocks = varargin{1}; % zeros restrictions on loadings
    ar_errors = varargin{2}; % logical: include ar(1) erros of the data
    threshold = varargin{3}; % threshold for likelihood function convergence
    varnames = (1:k)';
else
    blocks = varargin{1}; % zeros restrictions on loadings
    ar_errors = varargin{2}; % logical: include ar(1) erros of the data
    threshold = varargin{3}; % threshold for likelihood function convergence
    varnames = varargin{4};
    
end

% only used to print results at the end
varnames = array2table(varnames);
varnames.Properties.VariableNames = {'Series'};

% number of columns must equal the number of factors
if ~size(blocks,2) == m
    error('Size of blocks is not consistent with number of factors')
end

% number of rows must eaqual the numer of series k
if ~size(blocks,1) == k
    error('Size of blocks is not consistent with number of series')
end

if any(sum(blocks,1) == 1) && ar_errors
    error('With AR errors blocks must allow factors to load on more than one series')
end


%% Prepare data -----------------------------------------------------------
Mx = mean(X,'omitnan');
Wx = std(X,'omitnan');
% data including NaN values (i.e. not imputed)
xNaN = 10*(X-repmat(Mx,T,1))./repmat(Wx,T,1);  % Standardize series, ie scale()

% get the number of high frequnecy periods in every low frequency period
frq = set_frequencies(frq); 
%% Initial Conditions -----------------------------------------------------

[A_new, H_new, Q_new, R_new] = InitCond(xNaN,m,p,frq,isdiff,blocks);

% Initialize EM loop values
previous_loglik = -inf; % log likelihood can only get better
num_iter = 0; % value of the current iteration
converged = 0; % logical value: has the model converged?
max_iter = 1000; % maximum number of iterations

% Y for the estimation is WITH missing data
Y = xNaN'; %transpose for faster column-wise access

%% EM LOOP ----------------------------------------------------------------

%The model can be written as
%y = HJ*Z + e;
%Z = A*Z(-1) + v
% ~ means "not"
while (num_iter < max_iter) && ~converged % Loop until converges or max iter.
    
    H = H_new;
    R = R_new;
    A = A_new;
    Q = Q_new;
    
    [H_new, R_new, A_new, Q_new, loglik] = ...  % Applying EM algorithm
            EMstep(Y, A, H, Q, R, p, frq, isdiff, blocks);
    
    if num_iter > 2  % Checking convergence
        [converged, ~] = ...
            em_converged(loglik, previous_loglik, threshold, 1);
    end

    if (mod(num_iter,10) == 0) && (num_iter > 0)  % Print updates to command window
        disp(['Now running the ',num2str(num_iter),...
              'th iteration of max ', num2str(max_iter)]);
        disp(['  Loglik','   (% Change)'])
        disp([num2str(loglik),'   (', sprintf('%6.2f',100*((loglik-previous_loglik)/previous_loglik)) '%)'])
    end
    
    eA = abs(eig(A_new));
    
    if max(eA) >= 1
        disp('Estimated transition matrix non-stationary, breaking EM iterations')
        break
    end
    %LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;
    
    if converged
        H = H_new;
        R = R_new;
        A = A_new;
        Q = Q_new;
    end

end

if ar_errors
    
disp('Estimating model with AR errors')

[A_new, H_new, Q_new, R_new] = InitAR(Y,A,H,Q,R,frq,isdiff);

% Initialize EM loop values
previous_loglik = -inf;
num_iter = 0;
converged = 0;
    
while (num_iter < max_iter) && ~converged % Loop until converges or max iter.

    R = R_new;
    A = A_new;
    Q = Q_new;
    H = H_new;
    
    [H_new, R_new, A_new, Q_new, loglik] = ...  % Applying EM algorithm
            EMstepAR(Y, A, H, Q, R, p, frq, isdiff, blocks);
    
    if num_iter > 2  % Checking convergence
        [converged, ~] = ...
            em_converged(loglik, previous_loglik, threshold, 1);
    end

    if (mod(num_iter,10) == 0) && (num_iter > 0)  % Print updates to command window
        disp(['Now running the ',num2str(num_iter),...
              'th iteration of max ', num2str(max_iter)]);
        disp(['  Loglik','   (% Change)'])
        disp([num2str(loglik),'   (', sprintf('%6.2f',100*((loglik-previous_loglik)/previous_loglik)) '%)'])
    end
    
    eA = abs(eig(A_new));
    
    if max(eA) >= 1
        disp('Estimated transition matrix non-stationary, breaking EM iterations')
        break
    end
    %LL = [LL loglik];
    previous_loglik = loglik;
    num_iter =  num_iter + 1;
    
    if converged
        H = H_new;
        R = R_new;
        A = A_new;
        Q = Q_new;
    end

end

end

if(num_iter < max_iter)
    disp(['Successful: Convergence at ', num2str(num_iter), ' iterations'])
else
   disp('Stopped because maximum iterations reached')
end

% Final Run of filter/smoother using our estimated parameters
T = size(X_pred, 1);
HJ = get_HJ(H,frq,isdiff,p);
xpNaN = 10*(X_pred-repmat(Mx,T,1))./repmat(Wx,T,1);
if ar_errors
    HH = [HJ, eye(k)];
    [Zsmooth, Vsmooth, ~, LogLik, Update] = runKF(xpNaN', A, HH, Q, R);
else
    [Zsmooth, Vsmooth, ~, LogLik, Update] = runKF(xpNaN', A, HJ, Q, R);
end
Zsmooth = Zsmooth(:, 2:end)'; % Drop pre-sample values 
Vsmooth = Vsmooth(:, :, 2:end); % Drop pre-sample values 
var_Y = zeros(T,k);
if ar_errors
    sA = size(A,1);
    sa = sA - k; % size of A matrix excluding AR errors
    for t=1:T
        var_Y(t,:) = diag(HH*Vsmooth(:,:,t)*HH')';
    end
    y_common = Zsmooth(:,1:sa)*HH(:,1:sa)'; % common components (i.e. fit 
    % due to factors)
    y_ar = Zsmooth(:, sa+1:sA); % ar errors
    y_smooth = y_common + y_ar;  % Get smoothed Y
    Res.Y_common = repmat(Wx,T,1).*y_common/10 + repmat(Mx,T,1); 
    Res.Y_ar = repmat(Wx,T,1).*y_ar/10 + repmat(Mx,T,1);
else
    for t=1:T
        var_Y(t,:) = diag(HJ*Vsmooth(:,:,t)*HJ')' + R; % variance of predicted observations
    end
    y_smooth = Zsmooth * HJ';  % Get smoothed Y
end
% 1 s.d. confidence intervals ignoring parameter uncertainty
y_upper = y_smooth + sqrt(var_Y);
y_lower = y_smooth - sqrt(var_Y); 

%%  Loading the structure with the results --------------------------------
Res.y_smooth = y_smooth; %smoothed (fitted) values of inputs data
Res.Y_smooth = repmat(Wx,T,1).*y_smooth/10 + repmat(Mx,T,1);  % Unstandardized, smoothed values
Res.Y_upper = repmat(Wx,T,1).*y_upper/10 + repmat(Mx,T,1); % unstandardized upper bound
Res.Y_lower = repmat(Wx,T,1).*y_lower/10 + repmat(Mx,T,1); % unstandardized lower bound
Res.Z = Zsmooth; % factors
Res.H = H; % loadings
Res.HJ = HJ; % loadings in Kalman Filter format
Res.R = R; % Shocks in observation equatoin
Res.A = A; % Transition matrix 
Res.Q = Q; % Shocks in transition equation
Res.Mx = Mx; % Mean of data
Res.Wx = Wx; % scale parameter of data
Res.m = m; % number of factors
Res.p = p; % number of paremters
Res.forecast_loglikelihood = real(LogLik); % Log Likelihood for data used for prediction
Res.fitted_loglikelihood = real(loglik); % Log likelihood of data used for fitting the model
Res.Update = Update; % How much does each series contribute at each period

%% Display output
% Table with names and factor loadings

fprintf('\n\n\n');

try
disp('Table 1: Factor Loadings');
disp([varnames, array2table(Res.H)])  % Only select lag(0) terms
fprintf('\n\n\n');
catch
end

% Table with AR model on factors (factors with AR parameter and variance of residuals)

try
disp('Table 2: Observation Variance');
disp([varnames, array2table(Res.R')])  % Only select lag(0) terms
fprintf('\n\n\n');
catch
end

end
