# dfm-Matlab

Matlab functions to estimate and work with dynamic factor models (DFMs) developed in collaboration the the [BCC](https://www.bccprogramme.org)

## Description

This code is losely based on the open source code from "[Macroeconomic Nowcasting and Forecasting with Big Data](https://www.newyorkfed.org/research/staff_reports/sr830.html)" by Brandyn Bok, Daniele Caratelli, Domenico Giannone, Argia M. Sbordone, and Andrea Tambalotti, *Staff Reports 830*, Federal Reserve Bank of New York (prepared for Volume 10 of the *Annual Review of Economics*).

**Note:** This simplified code is written so that it should be easy to follow and implement for any standard mixed frequency data set. **These modifications were implemented by Seth Leonard in collaboration with the BCC and are not associated with the Federal Reserve Bank of New York or any of its staff.**

## Using this code
 
This code requires the [mutils](https://github.com/macroeconomicdata/mutils) toolbox. The function `dfm()` accepts mixed frequency data at the following frequencies:

 - 'd' daily
 - 'w' weekly
 - 'b' every other week
 - 'm' monthly
 - 'q' quarterly
 - 'y' yearly
 
See `dfm.m` and the examples in the scripts folder for further details. 

## Download instructions

Download the code as a ZIP file by clicking the green 'Clone or download' button and selecting 'Download ZIP'.

## File and folder description

* `functions/` : functions for estimating the model
* `scripts/load_process_DFM_free.m` : example script to estimate a dynamic factor model (DFM) for a panel of weekly and monthly data using free data from [macroeconomicdata.com](https://macroeconomicdata.com)
* `scripts/load_process_DFM_switzerland.m` : example script to estimate a dynamic factor model (DFM) for a panel of weekly and monthly data using Swiss data from [macroeconomicdata.com](https://macroeconomicdata.com)


## User Functions

- `Res = dfm(X,X_pred,m,p,frq,isdiff,blocks, threshold, ar_errors, varnames)` Main function for estimating dynamic factor models. The first six arguments are required; the remaining four are optional. 
- `S = KF(Y, A, HJ, Q, R)` Fast Kalman filtering adding each series sequentially (and thus avoiding matrix inversions).
- `S = KFilter(Y, A, HJ, Q, R)` Standard Kalman filtering entering data for each period as a vector. Results are identical to `KF()` above.
- `S = Ksmooth(A, S)` Kalman smoother. The input `S` is the output structure from either `KF()` or `KFilter()`.
- `[Zsmooth, Vsmooth, VVsmooth, loglik, Update] = runKF(Y, A, HJ, Q, R)` - Run the Kalman filter (`KF()`) and Kalman smoother.
- `[Y, X] = sim_dfm(T)` Simulate uniform frequency DFM data.
- `[Y, X] = sim_dfm_mixed_frq(T)` Simulated mixed frequency DFM data.

## Required software and versioning

MATLAB is required to run the code. The code was tested in MATLAB R2015b and later versions. Functionality with earlier versions of MATLAB is not guaranteed.
