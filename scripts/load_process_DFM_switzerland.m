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
macro_folder = 'C:\Users\seton\Dropbox\macroeconomicdata\macro';
addpath(macro_folder)

%% Load data for Switzerland
Xtable = readtable('swz_data_sa.csv');
% Check that the refference date is in fact formatted as a date:
class(Xtable.ref_date) % datetime
%% Load Library
lib = readtable('library.csv');
% Select only entries for the US:
lib = lib(lib.country == "switzerland", :);
% If you want to look at everything in there:
% lib.series_name
%% Select data
% For this example we'll use the following series:
selected = {'consumer price index cpi', 'import prices', 'currency',...
    'interbank rate', 'stock market', 'government bond 10y'};
% Get the data in which we're interested
cpi_data = Xtable(ismember(Xtable.series_name, selected), :);
% We'll also only use data since 1992:
cpi_data = cpi_data(cpi_data.ref_date >= '1995-01-01',:);

%% Process data
% This processes data using the information about the series in the library
% file. However, we'll need to tell the function whether to detrend data or
% standardize data. 
rm_trend = true;
standardize = true;
% Here we have monthly and daily data. To reduce the number of lags in the
% DFM, we will aggregate the daily data to weekly. 

lib(strcmp(lib.series_name, 'currency'),4) = {'weekly'};
lib(strcmp(lib.series_name, 'interbank rate'),4) = {'weekly'};
lib(strcmp(lib.series_name, 'stock market'),4) = {'weekly'};
lib(strcmp(lib.series_name, 'government bond 10y'),4) = {'weekly'};

% Xf excludes outliers
% Xp does not
% Estimate parameters excluding outliers, but you may want to make actual
% predictions using outliers if recent data was dropped.
% Note also that we are using the seasonally unadjusted series for CPI to
% provide an example of using seasonal adjustment in the process()
% function.

[Xf, Xp, Tnd, X_saf, X_level, X_scale, X_center] = ...
    process(cpi_data, lib, rm_trend, standardize);

% Because we are using weekly frequency for the data month dates
% have been indexed to the end of the week in which they fall

%% Plot processed data
% There's nothing better than eyeball econometrics. We'll plot the data
% without outliers, which is Xf; Xp includes outliers (f is for fit, p is
% for predict).
X = table2array(Xf);
% We'll have to fill in missing values to plot monthly data
X(:,1) = spline_fill_plot(X(:,1));
X(:,5) = spline_fill_plot(X(:,4));
figure(1)
plot(Xf.Time, X)
legend(Xf.Properties.VariableNames)

%% Nowcast the CPI: In Sample

% Print out the variables entering the model:
varnames = Xf.Properties.VariableNames'

% We'll construct a two factor model. Blocks allow us to put zero
% restrictions on the loadings of the model.
blocks = ones(6,2); 
% we enter a 0 where we want to impose a restriction. We'll (somewhat
% arbitrarily) make one block the CPI, currency, and stock market, and the
% other block CPI, governemnt bond yeilds, import prices, and the interbank
% rate. Not that blocks are optional; the default is no restrictions.
blocks([3,4,5],1) = 0; % currency and stock market block
blocks([2,6],2) = 0; %  import prices, and interbank block
m = 2; % this model has 2 factors
p = 4; % we'll include 4 lags in the transition equation
% Note the order of these variables MUST be the same as the column order of
% the data Xf and Xp.
frq = lib{ismember(lib.series_name, selected), 4}; % frequency of each
% input (not case sensitive). dfm() only uses the first letter, so that you
% can enter, for example {'m', 'w', 'q'} for month, week, quarter.
isdiff = lib{ismember(lib.series_name, selected), 9};
threshold = 1e-5; % threshold for EM convergence
ar_errors = true; % allow for AR(1) error terms in variables
X = table2array(Xf);
X_pred = table2array(Xp);

% Estimate the DFM
Res = dfm(X,X_pred,m,p,frq,isdiff, blocks, ar_errors, threshold, varnames);

%% Plot Fitted Values
% Fill missing weekly values in the true data for plotting
true_vals = spline_fill_plot(X_pred(:,1));
figure(1) % Plot the results with one standard deviation confidence intervals
h = plot([true_vals, Res.Y_smooth(:,1), Res.Y_upper(:,1), Res.Y_lower(:,1)]);
set(h, {'color'}, {'black'; 'red'; 'blue'; 'blue'})
set(h, {'lineStyle'}, {'-'; '-'; '--'; '--'})
legend({'Detrended CPI Growth', 'Fitted Value'})

%% Decomposing the CPI 
% We can decompose the CPI into three parts: a common component (i.e. due 
% to the factors of the model), an AR(1) component, and a residual, which
% we will call y_iid below.
y_ar = Res.Y_ar(:,1); % AR(1) component
y_common = Res.Y_common(:,1); % common component
y_iid = X_pred(:,1) - y_ar - y_common; % redidual
y_smooth = Res.y_smooth(:,1); % This is the AR(1) plus the common component

% In this case volatility in the CPI (i.e. log differences in the series)
% is due almost entierely to the AR(1) and IID components
figure(2)
bar([y_ar, y_common, y_iid], 'stacked')
legend({'AR errors', 'Common', 'IID'})

%% Converting results back to levels
Y_pred = Res.Y_smooth; % predictions for all variables in model

% To use unprocess data must be a timetable and column names must match
% names in the lib file
Xtimetable = array2timetable(Y_pred, 'RowTimes', Xp.Time);
Xtimetable.Properties.VariableNames = Xp.Properties.VariableNames;

% Convert everything back to levels
% X_sa is seasonally adjusted (not adding back in seasonal components)
% X_nsa is not seasonally adjusted; seasonal components are added back in
% If we did not seasonally adjust the data, these two will be the same
[X_sa, X_nsa] = unprocess(Xtimetable, lib, X_level, Tnd, X_saf, X_scale, X_center);

% We can plot predicted against actual values as follows:
X_true = Xtable(strcmp(Xtable.series_name, 'consumer price index cpi'), :);
X_true = X_true(X_true.ref_date >= '1995-01-01',:); % we used data from 1995
X_true.ref_date = dateshift(X_true.ref_date, 'end', 'week'); % because our
% model is weekly, we will have to index dates to the end of the week.
X_true = timetable(X_true.ref_date, X_true.value); % convert data to timetable

cpi_levels = innerjoin(X_true, X_sa(:,1)); % The fact that these are timetables 
cpi_levels = innerjoin(cpi_levels, X_nsa(:,1));
% insures that dates will line up correctly
figure(3)
plot(cpi_levels.Time, cpi_levels{:,:})
legend({'True Values', 'Fitted SA Values', 'Fitted NSA Values'})




