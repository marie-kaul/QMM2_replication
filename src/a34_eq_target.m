function res = a34_eq_target(vars_A33, vars_A32, theta_o, psi, targets, params)

%========================================================================%
% A.3.4 Solving for psi and theta_o
%
% Inputs:
%   x(1) = psi        : fraction of investors among buyers
%   x(2) = theta_o    : ownership-market tightness
%
%   vars   : struct with variables from A.3.1–A.3.3
%   params : struct with parameters
%
% Output:
%   res(1) : residual of equation (A.35)
%   res(2) : residual of equation (A.74)
%========================================================================%

%% ----------------------------------------------------------
% 0a. Unpack parameters
% ----------------------------------------------------------

r       = params.r;
rho     = params.rho;

Fh      = params.Fh;
Fl      = params.Fl;

omega_h = params.omega_h;
omega_l = params.omega_l;

%% ----------------------------------------------------------
% 0b. Unpack variables from previous sections
% ----------------------------------------------------------

% Market tightness & meeting rates
theta_l = vars_A33.theta_l;
vo      = vars_A32.vo;
vl      = vars_A33.v_l;

% Surpluses
Sigma_h = vars_A32.Sigma_h;
Sigma_l = vars_A33.Sigma_l;

% Thresholds
Z       = vars_A33.Z;

% Housing stocks on the market
u_o     = vars_A32.u_o;
u_l     = vars_A33.u_l;

n = vars_A33.n;

%% ----------------------------------------------------------
% 1. Equation (A.35): population / flow balance
% ----------------------------------------------------------

% Buyers in ownership and rental markets
buyers_o = vo * theta_o * u_o;
buyers_l = vl * theta_l * u_l;

% Fraction psi of buyers are investors
lhs_A35 = psi * buyers_o + (1 - psi) * buyers_l;

% Total population
rhs_A35 = n;

res_A35 = lhs_A35 - rhs_A35;

%% ----------------------------------------------------------
% 2. Equation (A.74): marginal buyer indifference
% ----------------------------------------------------------
lhs_A74 = (1 - omega_h) * vo * Sigma_h ...
        - (1 - omega_l) * vl * Sigma_l;

rhs_A74 = (r + rho) * Z + Fh - Fl;

res_A74 = lhs_A74 - rhs_A74;

%% ----------------------------------------------------------
% 6. Collect residuals
% ----------------------------------------------------------
res = [res_A35; res_A74];

end
