% First we will clear up our workspace; make sure you have saved everything
% you need!
close all
clear 
%% Add paths for functions and data
% These directories contains functions we will use from 
% https://github.com/macroeconomicdata
addpath 'mutils/functions' 
addpath 'dfm-matlab/functions'
addpath(genpath('mutils/toolbox')) % this is where the X-13 toolbox is located
% Location of macroeconomicdata.com macro dropbox folder
macro_folder = 'C:\Users\seton\Dropbox\sample_data'; % Get the data
addpath(macro_folder)

%% Load free data
Xtable = readtable('usa_free.csv');
% Check that the refference date is in fact formatted as a date:
class(Xtable.ref_date) % datetime
%% Load Library
lib = readtable('library_free.csv');
% If you want to look at everything in there:
% lib.series_name
%% Select data
% For this example we'll model personal consumption expenditures using
% the following series:
selected = {'personal consumption expenditures',...
    'personal consumption expenditures durable goods',...
    'personal consumption expenditures services' ,...
    'mortgage spread 30y 15y', 'initial jobless claims',...
    'gasoline prices sa', 'bank lending rate', 't bill spread 10y 3m'};
% Get the data in which we're interested
pce_data = Xtable(ismember(Xtable.series_name, selected), :);
% We'll also only use data since 1992:
pce_data = pce_data(pce_data.ref_date >= '1991-01-01',:);

%% Select data
% This processes data using the information about the series in the library
% file. However, we'll need to tell the function whether to detrend data or
% standardize data. 
rm_trend = true;
standardize = true;
% Here we have monthly daily, and weekly data. To reduce the number of lags
% in the DFM, we will aggregate the daily data to weekly. To do so we will
% set the 'frequency' in lib for the daily series to 'weekly'

lib(strcmp(lib.series_name, 't bill spread 10y 3m'),4) = {'weekly'};
lib(strcmp(lib.series_name, 'bank lending rate'),4) = {'weekly'};

% To take a look at the library entries we will use:
lib(ismember(lib.series_name, selected),:)

%% Process data
% Xf excludes outliers
% Xp does not
% Estimate parameters excluding outliers, but you may want to make actual
% predictions using outliers if recent data was dropped.
% Note also that we are using the seasonally unadjusted series for CPI to
% provide an example of using seasonal adjustment in the process()
% function.

[Xf, Xp, Tnd, X_saf, X_level, X_scale, X_center] = ...
    process(pce_data, lib, rm_trend, standardize);

% Because we are using weekly frequency for the data month dates
% have been indexed to the end of the week in which they fall

%% Plot processed data
% There's nothing better than eyeball econometrics. We'll plot the data
% without outliers, which is Xf; Xp includes outliers (f is for fit, p is
% for predict).
X = table2array(Xf);
% We'll have to fill in missing values to plot monthly data
X(:,2) = spline_fill_plot(X(:,2));
X(:,3) = spline_fill_plot(X(:,3));
X(:,6) = spline_fill_plot(X(:,6));
X(:,7) = spline_fill_plot(X(:,7));
X(:,8) = spline_fill_plot(X(:,8));
figure(1)
plot(Xf.Time, X)
legend(Xf.Properties.VariableNames)

%% Nowcast the PCE: In Sample

% Print out the variables entering the model:
varnames = Xf.Properties.VariableNames';

% We'll construct a two factor model. Blocks allow us to put zero
% restrictions on the loadings of the model.
blocks = ones(8,2); 
% we enter a 0 where we want to impose a restriction. We'll (somewhat
% arbitrarily) make one block the CPI, currency, and stock market, and the
% other block CPI, governemnt bond yeilds, import prices, and the interbank
% rate. Not that blocks are optional; the default is no restrictions. In
% that case the model will be under-identified and estimation with AR error
% terms is not recommended.
blocks([1,2,4,8],1) = 0; % real block
blocks(3,2) = 0; %  prices block
m = 2; % this model has 3 factors
p = 4; % we'll include 4 lags in the transition equation
% Note the order of these variables MUST be the same as the column order of
% the data Xf and Xp.
frq = lib{ismember(lib.series_name, selected), 4}; % frequency of each
% input (not case sensitive). dfm() only uses the first letter, so that you
% can enter, for example {'m', 'w', 'q'} for month, week, quarter.
isdiff = lib{ismember(lib.series_name, selected), 9};
threshold = 1e-5; % threshold for EM convergence
ar_errors = false; % do not allow for AR(1) error terms in variables (for 
% this model including AR(1) errors causes the iterations to fail to
% converge due to weak identification). 
X = table2array(Xf);
X_pred = table2array(Xp);

% Estimate the DFM
Res = dfm(X,X_pred,m,p,frq,isdiff, blocks, ar_errors, threshold, varnames);

% The model overfitts 't bill spread 10y 3m', so that the prices
% factor is simply this observed series.
%% Plot Fitted Values
% Fill missing weekly values in the true data for plotting
id = 5;
true_vals = spline_fill_plot(X_pred(:,id));
figure(1) % Plot the results with one standard deviation confidence intervals
h = plot([true_vals, Res.Y_smooth(:,id), Res.Y_upper(:,id), Res.Y_lower(:,id)]);
set(h, {'color'}, {'black'; 'red'; 'blue'; 'blue'})
set(h, {'lineStyle'}, {'-'; '-'; '--'; '--'})
legend({'Detrended PCE Growth', 'Fitted Value'})

%% Decomposing factor contributions to PCE

z = Res.Z; % factors
h = Res.HJ(id,:); % loadings on PCE
HH = kron(ones(size(Xf,1),1),h); % expanding loadings for element-wise
% multiplication
F = Res.Wx(id)*HH.*z + Res.Mx(id); % Contribution of each factor to PCE 
% adjusted for scaling/centering (Wx is the scale parameter and Mx the 
% centering parameter
f1 = sum(F(:,1:m:end),2); % first factor (summing all lags)
f2 = sum(F(:,2:m:end),2); % second factor (summing all lags)

% Check that everything adds up!
max(f1 + f2 - Res.Y_smooth(:,id)) % zero

figure(4)
hold on
bar([f1,f2], 'stacked');
plot(Res.Y_smooth(:,id), 'color', 'red');
hold off
legend({'real', 'prices'})

% The large magnitudes of each factor's contribution in oposite directions
% is indicative of overfitting.

%% Converting results back to levels
% To use unprocess data must be a timetable and column names must match
% names in the lib file
Xtimetable = array2timetable(Res.Y_smooth, 'RowTimes', Xp.Time);
Xtimetable.Properties.VariableNames = Xp.Properties.VariableNames;

% Convert everything back to levels
% X_sa is seasonally adjusted (not adding back in seasonal components)
% X_nsa is not seasonally adjusted; seasonal components are added back in
% If we did not seasonally adjust the data, these two will be the same
[X_sa, X_nsa] = unprocess(Xtimetable, lib, X_level, Tnd, X_saf, X_scale, X_center);

% We can plot predicted against actual values as follows:
X_true = Xtable(strcmp(Xtable.series_name, 'personal consumption expenditures'), :);
X_true = X_true(X_true.ref_date >= '1991-01-01',:); % we used data from 1995
X_true.ref_date = dateshift(X_true.ref_date, 'end', 'week'); % because our
% model is weekly, we will have to index dates to the end of the week.
X_true = timetable(X_true.ref_date, X_true.value); % convert data to timetable

pce_levels = innerjoin(X_true, X_sa(:,id)); % The fact that these are timetables 
% insures that dates will line up correctly
figure(3)
plot(pce_levels.Time, pce_levels{:,:})
legend({'True Values', 'Fitted Values'})




