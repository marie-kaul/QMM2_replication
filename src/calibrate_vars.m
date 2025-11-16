%========================================================================%
% This file calibrates the exogenous variables.
% Input: nothing.
% Output: data/calibration_targets.mat
%========================================================================%

clear; clc;

%% ----------------------------------------------------------
% A. Directly imposed targets
%% ----------------------------------------------------------

n   = 1;          % Equal numbers of households and properties
chi = 1;          % Speed of population adjustment
Bl  = 0;          % No incentive for households to leave

omega_ratio = 1;  % ωo/ηo = ωl/ηl

Fh_vo_over_Y = 0.5;      % (Fh/vo)/(Y/365)
Lambda_l_over_Lambda_h = 0.5;  % Λl / Λh
Fl_vl_over_Fh_vo = 0.5;        % (Fl/vl)/(Fh/vo)
Fi_over_Fh = 1;                % Fi / Fh

%% ----------------------------------------------------------
% B. Empirical targets
%% ----------------------------------------------------------

h = 0.54;              % Homeownership rate
i_share = 0.054;       % Buy-to-rent share
Pi_over_R  = 14.5;     % Price-rent ratio
Pi_over_Ph = 0.99;     % Investor/home-buyer price ratio

phi_first = 0.40;      % Fraction of first-time buyers
Xi_age = 8.3;          % Age difference

rg   = 0.0186;         % Risk-free real interest rate
rkbar = 0.0493;        % Average real mortgage interest rate
rz   = 0.0643;         % Marginal real mortgage rate

l_LTV = 0.80;          % Loan-to-value ratio
Tk = 25;               % Mortgage term

Ch_over_Ph = 0.00;     % Buyer non-tax transaction cost
Ci_over_Pi = 0.00;     % Investor non-tax buyer cost (same)

D_over_P = 0.026;      % Maintenance cost
Dl_over_R = 0.08;      % Landlord extra maintenance
Co_over_P = 0.045;     % Seller transaction cost
Cl_over_R = 0.083;     % Landlord transaction cost
A_over_Cl = 0.00;      % Tenancy agreement fee

Tso = 0.161;           % Sellers’ average time on market
Tbh = 0.206;           % Home-buyers’ average time on market
Tsl = 0.066;           % Landlords’ time-to-lease

Lambda_h = 20.6;       % Viewings per home-buyer

Tmh = 9.25;            % Time between moves (owners)
Tml = 3.04;            % Time between moves (tenants)

Ph_over_Y = 5.6;       % House price / income ratio
P_avg = 402e3;         % Average transaction price ($402k)

tau_h = 0.015;         % LTT rate
tau_i = 0.015;         % Same for investors

%% ----------------------------------------------------------
% B. Matched response to new LTT
%% ----------------------------------------------------------

beta_mh = -0.13;       % Change in log moving rate of owners

%% ----------------------------------------------------------
% Package into a struct for convenience
%% ----------------------------------------------------------

targets = struct( ...
    'n', n, 'chi', chi, 'Bl', Bl, 'omega_ratio', omega_ratio, ...
    'Fh_vo_over_Y', Fh_vo_over_Y, ...
    'Lambda_l_over_Lambda_h', Lambda_l_over_Lambda_h, ...
    'Fl_vl_over_Fh_vo', Fl_vl_over_Fh_vo, ...
    'Fi_over_Fh', Fi_over_Fh, ...
    'h', h, 'i_share', i_share, ...
    'Pi_over_R', Pi_over_R, 'Pi_over_Ph', Pi_over_Ph, ...
    'phi_first', phi_first, 'Xi_age', Xi_age, ...
    'rg', rg, 'rkbar', rkbar, 'rz', rz, ...
    'l_LTV', l_LTV, 'Tk', Tk, ...
    'Ch_over_Ph', Ch_over_Ph, 'Ci_over_Pi', Ci_over_Pi, ...
    'D_over_P', D_over_P, 'Dl_over_R', Dl_over_R, ...
    'Co_over_P', Co_over_P, 'Cl_over_R', Cl_over_R, ...
    'A_over_Cl', A_over_Cl, ...
    'Tso', Tso, 'Tbh', Tbh, 'Tsl', Tsl, ...
    'Lambda_h', Lambda_h, ...
    'Tmh', Tmh, 'Tml', Tml, ...
    'Ph_over_Y', Ph_over_Y, 'P_avg', P_avg, ...
    'tau_h', tau_h, 'tau_i', tau_i, ...
    'beta_mh', beta_mh ...
);

%% ----------------------------------------------------------
% Save calibration targets to a .mat file in ../data/
%% ----------------------------------------------------------

% Path to target folder (one level up from src, into data/)
save_folder = fullfile('..', 'data');

% Create folder if it does not exist
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Full path to file
save_path = fullfile(save_folder, 'calibration_targets.mat');

% Save all relevant variables or the struct 'targets'
save(save_path, 'targets');

fprintf('Calibration data saved to: %s\n', save_path);


