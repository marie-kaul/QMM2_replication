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

load('../data/calibration_targets.mat');

%% ----------------------------------------------------------
% A.3.1 Ownership-market thresholds
%% ----------------------------------------------------------

% Guess psi and theta_o for now. Later: outer loop over different values.
psi = 0.003;
theta_o = 1;


