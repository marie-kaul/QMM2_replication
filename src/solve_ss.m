%========================================================================%
% This file solves for steady state.
% Input: data/calibration_targets.mat.
% Output: 
% 
% From Section A.3 of the paper: 
% The solution method is based on a numerical search over the fraction ψ of 
% investors among buyers and ownership-market tightness θo that satisfy two 
% equations representing equilibrium in the ownership and rental markets. 
% Within this search, given a (ψ,θo), the ownership-market thresholds (yh, xh) 
% and rental-market and credit-cost thresholds (yl ,Z) are found by solving 
% two equations numerically.
%========================================================================%

clear; clc;

load('../data/calibration_data.mat');

%% ----------------------------------------------------------
% Initial variable guesses
%% ----------------------------------------------------------

psi = 0.003;
theta_o = 1;

%% ----------------------------------------------------------
% A.3.1 Ownership-market thresholds
%% ----------------------------------------------------------

vars_A31 = a31_ownership_mkt(psi, theta_o, targets, params);

yh = vars_A31.yh;
xh = vars_A31.xh;

if params.delta_h*yh < xh 
    disp('A.3.1: δhyh < xh is satisfied.')
else
    error('A.3.1: δhyh < xh is not satisfied.')
end

%% ----------------------------------------------------------
% A.3.2 Ownership-market thresholds
%% ----------------------------------------------------------

vars_A32 = a32_ownership_other(yh,xh, psi, theta_o, targets, params);

%% ----------------------------------------------------------
% A.3.3 The rental-market and credit-cost thresholds, 
% and the city population
%% ----------------------------------------------------------

vars_A33 = a33_rental_credit_population(vars_A31, vars_A32, theta_o, ...
    psi, targets, params);

