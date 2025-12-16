%========================================================================%
%
% Constructs calibration targets (Table 3) and calibrated parameters
% (Table 4).
%
% Output: ../data/calibration_data.mat
%         contains structs:
%           - targets : all calibration targets from Table 3
%           - params  : all calibrated parameter values from Table 4
%
% NOTE ON UNITS (OPTION B — OFFICIAL REPLICATION CONVENTION):
%   - Time units: years
%   - Monetary units: THOUSANDS of 2007 Canadian dollars
%       * "$10.5k" in the paper is stored as 10.5
%       * "$402k" is stored as 402
%   - Credit-cost distribution:
%       log(K) ~ N(mu, sigma^2), where K is measured in *thousands* of dollars
%       → mu = 5.05 is used directly (NO log(1000) adjustment)
%
% IMPORTANT:
%   All equations must be interpreted in thousands-of-dollars units.
%========================================================================%

clear; clc;

pwd
which calibration_data.mat -all


%% ----------------------------------------------------------
% A. Targets for pre-tax-change steady state (Table 3)
%% ----------------------------------------------------------

targets = struct();

% -------------------------------
% Directly imposed targets
% -------------------------------
targets.chi                = 1;
targets.Bl                 = 0;
targets.omega_ratio        = 1;

targets.Fh_vo_over_Y_day        = 0.5;
targets.Lambda_l_over_Lambda_h  = 0.5;
targets.Fl_vl_over_Fh_vo        = 0.5;
targets.Fi_over_Fh              = 1;

% -------------------------------
% Empirical targets
% -------------------------------
targets.h                  = 0.54;
targets.i_share            = 0.054;

targets.Pi_over_R          = 14.5;
targets.Pi_over_Ph         = 0.99;

targets.phi_first          = 0.40;
targets.Xi_age             = 8.3;

targets.rg                 = 0.0186;
targets.rkbar              = 0.0493;
targets.rz                 = 0.0643;

targets.l_LTV              = 0.80;
targets.Tk                 = 25;

targets.Ch_over_Ph         = 0.0;
targets.Ci_over_Pi         = 0.0;

targets.D_over_P           = 0.026;
targets.Dl_over_R          = 0.08;

targets.Co_over_P          = 0.045;
targets.Cl_over_R          = 0.083;
targets.A_over_Cl          = 0.0;

targets.Tso                = 0.161;
targets.Tbh                = 0.206;
targets.Tsl                = 0.066;

targets.Lambda_h           = 20.6;

targets.Tmh                = 9.25;
targets.Tml                = 3.04;

targets.Ph_over_Y          = 5.6;
targets.P_avg              = 402;      % $402k → 402 (thousands)

targets.tau_h              = 0.015;
targets.tau_i              = 0.015;

% -------------------------------
% Matched response to new LTT
% -------------------------------
targets.beta_mh            = -0.13;


%% ----------------------------------------------------------
% B. Calibrated parameters (Table 4)
%     ALL MONETARY VALUES IN THOUSANDS OF DOLLARS
%% ----------------------------------------------------------

params = struct();

% Discounting and exit
params.r       = 0.0328;
params.rho     = 0.0427;
params.rho_l   = 0.00700;

% Maintenance costs
params.D       = 10.5;
params.Dl      = 2.20;

% Minimum match quality
params.zeta_h  = 33.6;
params.zeta_l  = 24.6;

% Match quality distributions
params.lambda_h = 33.1;
params.lambda_l = 36.2;

% Match quality shocks
params.alpha_h  = 0.0793;
params.alpha_l  = 0.279;
params.delta_h  = 0.855;

% Credit cost distribution (log of *thousands* of dollars)
params.xi       = 0.0828;
params.mu       = 5.05;      % <-- DO NOT rescale
params.sigma    = 0.674;

% Transaction costs
params.Ci       = 0.0;
params.Ch       = 0.0;

params.Co       = 18.1;
params.Cl       = 2.28;
params.Cw       = 0.709;

% Flow search costs
params.Fh       = 9.84;
params.Fi       = 9.84;
params.Fl       = 9.72;

% Entry cost to the city
params.E        = 21.1;

% Population adjustment
params.chi      = targets.chi;

% Meeting-function productivities
params.nu_o     = 110;
params.nu_l     = 165;

% Meeting-function elasticities
params.eta_o    = 0.490;
params.eta_l    = 0.764;

% Bargaining powers
params.omega_h  = 0.490;
params.omega_i  = 0.265;
params.omega_l  = 0.764;

% Convenience: store tax rates
params.tau_h    = targets.tau_h;
params.tau_i    = targets.tau_i;


%% ----------------------------------------------------------
% C. Save to ../data/ (ROBUST TO WORKING DIRECTORY)
%% ----------------------------------------------------------

this_file = mfilename('fullpath');
this_dir  = fileparts(this_file);

save_folder = fullfile(this_dir, '..', 'data');

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

save_path = fullfile(save_folder, 'calibration_data.mat');
save(save_path, 'targets', 'params');

fprintf('Saved calibration data to:\n%s\n', save_path);
