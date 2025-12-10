%========================================================================%
%
% Constructs calibration targets (Table 3) and calibrated parameters
% (Table 4). 
%
% Output: ../data/calibration_data.mat
%         contains structs:
%           - targets : all calibration targets from Table 3
%           - params  : all calibrated parameter values from Table 4
%========================================================================%

clear; clc;

%% ----------------------------------------------------------
% A. Targets for pre-tax-change steady state (Table 3)
%% ----------------------------------------------------------

% -------------------------------
% Directly imposed targets
% -------------------------------
targets.n                  = 1;        % Equal numbers of households and properties in the city
targets.chi                = 1;        % Speed of adjustment of the city population
targets.Bl                 = 0;        % No incentive for households to choose to leave the city

% Bargaining powers equal to meeting-function elasticities: ωo/ηo = ωl/ηl = 1
targets.omega_ratio        = 1;

% Cost per viewing for home-buyers relative to daily income (Fh/vo)/(Y/365)
targets.Fh_vo_over_Y_day   = 0.5;

% Viewings per renter relative to viewings per home-buyer Λl/Λh
targets.Lambda_l_over_Lambda_h = 0.5;

% Cost of a rental viewing relative to a home-buyer viewing (Fl/vl)/(Fh/vo)
targets.Fl_vl_over_Fh_vo   = 0.5;

% Flow search costs of investors relative to home-buyers Fi/Fh
targets.Fi_over_Fh         = 1;

% -------------------------------
% Empirical targets
% -------------------------------

targets.h                  = 0.54;     % Homeownership rate (54%)
targets.i_share            = 0.054;    % Buy-to-rent as a share of all transactions (5.4%)

targets.Pi_over_R          = 14.5;     % Average price-rent ratio for the same properties
targets.Pi_over_Ph         = 0.99;     % Average price paid by investors relative to home-buyers (99%)

targets.phi_first          = 0.40;     % Fraction of first-time buyers among all home-buyers (40%)
targets.Xi_age             = 8.3;      % Difference of average ages of owner-occupiers and tenants

targets.rg                 = 0.0186;   % Risk-free real interest rate (1.86)
targets.rkbar              = 0.0493;   % Average real mortgage interest rate r̄k (4.93%)
targets.rz                 = 0.0643;   % Real mortgage rate faced by marginal home-buyer (6.43%)

targets.l_LTV              = 0.80;     % Initial loan-to-value ratio of first-time buyers (80%)
targets.Tk                 = 25;       % Mortgage term (years)

% Non-tax transaction costs of buyers as a fraction of price Ch/Ph = Ci/Pi (0%)
targets.Ch_over_Ph         = 0.0;
targets.Ci_over_Pi         = 0.0;

targets.D_over_P           = 0.026;    % Property maintenance costs as a fraction of price (2.6%)
targets.Dl_over_R          = 0.08;     % Landlords’ extra maintenance/management costs relative to rent (8%)

targets.Co_over_P          = 0.045;    % Sellers’ transaction costs as a fraction of price (4.5%)
targets.Cl_over_R          = 0.083;    % Landlords’ transaction costs as a fraction of rent (8.3%)
targets.A_over_Cl          = 0.0;      % New tenancy agreement fee as a fraction of landlord costs (0%)

targets.Tso                = 0.161;    % Sellers’ average time on the market (years)
targets.Tbh                = 0.206;    % Home-buyers’ average time on the market (years)
targets.Tsl                = 0.066;    % Landlords’ average time on the rental market (years)

targets.Lambda_h           = 20.6;     % Average viewings per home-buyer

targets.Tmh                = 9.25;     % Average time between moves for owner-occupiers (years)
targets.Tml                = 3.04;     % Average time between moves for tenants (years)

targets.Ph_over_Y          = 5.6;      % Ratio of house prices to income Ph/Y
targets.P_avg              = 402e3;    % Average transaction price of a property ($402k)

targets.tau_h              = 0.015;    % Effective land transfer tax rate for all buyers (1.5%)
targets.tau_i              = 0.015;    % Same rate for investors (for completeness)

% -------------------------------
% B. Matched response to new LTT
% -------------------------------

targets.beta_mh            = -0.13;    % Change in log moving rate of owner-occupiers


%% ----------------------------------------------------------
% B. Calibrated parameters (Table 4)
%% ----------------------------------------------------------

params = struct();

% Time units are years; monetary units are 2007 Canadian dollars.
% Percentages in the table are stored as decimals here.

% Discounting and exit
params.r       = 0.0328;     % Discount rate for future housing-market payoffs (3.28%)
params.rho     = 0.0427;     % Households’ exit rate from the city (4.27%)
params.rho_l   = 0.00700;    % Investors’ exit rate (0.700%)

% Costs and match-utility levels (k = 1e3)
params.D       = 10.5e3;     % Property maintenance cost ($10.5k)
params.Dl      = 2.20e3;     % Landlords’ extra maintenance/management costs ($2.20k)

params.zeta_h  = 33.6e3;     % Minimum new match quality in the ownership market ($33.6k)
params.zeta_l  = 24.6e3;     % Minimum new match quality in the rental market ($24.6k)

% Match quality distributions (Pareto-type)
params.lambda_h = 33.1;      % Shape parameter of home-buyer new match quality distribution
params.lambda_l = 36.2;      % Shape parameter of tenant new match quality distribution

% Match quality shocks
params.alpha_h  = 0.0793;    % Arrival rate of match quality shocks (ownership market) (7.93%)
params.alpha_l  = 0.279;     % Arrival rate of match quality shocks (rental market) (27.9%)
params.delta_h  = 0.855;     % Size of match quality shocks in the ownership market

% Credit cost distribution
params.xi       = 0.0828;    % Fraction of tenants drawing a new credit cost after a shock (8.28%)
params.mu       = 5.05;      % Mean of distribution of credit costs
params.sigma    = 0.674;     % Std dev of distribution of credit costs

% Transaction costs (levels)
params.Ci       = 0.0;       % Buyers’ non-tax transaction costs (Ci = Ch = 0)
params.Ch       = 0.0;       % (store explicitly as well)

params.Co       = 18.1e3;    % Sellers’ transaction costs ($18.1k)
params.Cl       = 2.28e3;    % Landlords’ transaction costs ($2.28k)
params.Cw       = 0.709e3;   % Tenants’ transaction costs excl. tenancy-agreement fee ($0.709k)

% Flow search costs
params.Fh       = 9.84e3;    % Flow search costs of home-buyers ($9.84k)
params.Fi       = 9.84e3;    % Flow search costs of investors (Fi = Fh)
params.Fl       = 9.72e3;    % Flow search costs of prospective tenants in rental market ($9.72k)

% Entry cost to the city
params.E        = 21.1e3;    % Entry cost of moving to the city ($21.1k)

% Population adjustment
params.chi      = 1;         % Speed of adjustment of the city population (matches Table 3)

% Meeting-function productivities
params.nu_o     = 110;       % Ownership-market meeting productivity parameter
params.nu_l     = 165;       % Rental-market meeting productivity parameter

% Meeting-function elasticities (Hosios-type condition)
params.eta_o    = 0.490;     % Elasticity wrt sellers in ownership market
params.eta_l    = 0.764;     % Elasticity wrt landlords in rental market

% Bargaining powers
params.omega_h  = 0.490;     % Seller meeting a home-buyer
params.omega_i  = 0.265;     % Seller meeting an investor
params.omega_l  = 0.764;     % Landlord meeting a prospective tenant

% (Optionally) store tax rates here too, for convenience
params.tau_h    = targets.tau_h;
params.tau_i    = targets.tau_i;


%% ----------------------------------------------------------
% C. Save to ../data/
%% ----------------------------------------------------------

save_folder = fullfile('..', 'data');

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Main file with both targets and params
save_path = fullfile(save_folder, 'calibration_data.mat');
save(save_path, 'targets', 'params');
