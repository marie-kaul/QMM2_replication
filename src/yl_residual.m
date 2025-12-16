function vars = yl_residual(yl, vars_A31, vars_A32, ...
    theta_o, psi, targets, params)

%% ----------------------------------------------------------
% 0. Unpack variables from A.3.1–A.3.2
%% ----------------------------------------------------------

Sigma_i = vars_A31.Sigma_i;

Sigma_h = vars_A32.Sigma_h;
pi_h    = vars_A32.pi_h;
m_h     = vars_A32.m_h;
q_h     = vars_A32.q_h;
u_o     = vars_A32.u_o;
vo      = vars_A32.vo;
G       = vars_A32.G;

%% Parameters

r      = params.r;
rho    = params.rho;
rho_l  = params.rho_l;
sigma  = params.sigma;

zeta_l   = params.zeta_l;
lambda_l = params.lambda_l;

alpha_l = params.alpha_l;
ml      = alpha_l + rho_l;

Dl = params.Dl;
Cl = params.Cl;
Cw = params.Cw;
Fl = params.Fl;

D  = params.D;
Co = params.Co;
Ci = params.Ci;

nu_l  = params.nu_l;
eta_l = params.eta_l;

omega_l = params.omega_l;
omega_h = params.omega_h;
omega_i = params.omega_i;

xi   = params.xi;
mu   = params.mu;

E     = params.E;

% --- tax-adjusted bargaining powers ---
tau_i = targets.tau_i;
tau_h = targets.tau_h;

omega_h_star = omega_h / (1 + tau_h * (1 - omega_h));
omega_i_star = omega_i / (1 + tau_i * (1 - omega_i));


%% ----------------------------------------------------------
% Section: Rental-Market variables conditional on the transactions
% threshold 
%% ----------------------------------------------------------

% 1) Transaction probability pi_l (from (40))
pi_l = (zeta_l / yl)^lambda_l;

% 2) Surplus Sigma_l = pi_l*yl / ((lambda_l-1)(r+rho+ml))
Sigma_l = (pi_l * yl) / ((lambda_l - 1) * (r + rho + ml));

% 3) Letting rate s_l implied by yl via (A.64)
TERM1 = (1 + tau_i * (1 + rho_l/r)) ...
        * theta_o * vo ...
        * ((1 - psi) * omega_h_star * Sigma_h + psi * omega_i_star * Sigma_i);

TERM2 = (r + rho_l) * ((1 + tau_i) * Co + Ci + (1 + tau_i * omega_i_star) * Sigma_i);

TERM3 = - tau_i * (1 + rho_l/r) * D;

RHS_A63 = TERM1 + TERM2 + TERM3;

s_l = ((lambda_l - 1) * (r + rho + ml)) / (omega_l * yl) * RHS_A63;


% 4) Rental market tightness theta_l from (A.65)
theta_l = (s_l / (nu_l * pi_l))^(1 / (1 - eta_l));

% 5) Meeting rate v_l from (43)
v_l = nu_l * theta_l^(-eta_l);

%% ----------------------------------------------------------
% Section: Solving for the credit-cost threshold and city population
%% ----------------------------------------------------------

% Rental stocks (u_l, q_l) from steady-state flow conditions
den = 1 + s_l / (ml + rho);

u_l = (1 - q_h - u_o) / den;

q_l = (s_l / (ml + rho)) * u_l;

% Rental market objects
vars.theta_l = theta_l;
vars.v_l     = v_l;
vars.Sigma_l = Sigma_l;
vars.u_l     = u_l;
vars.q_l     = q_l;
vars.s_l     = s_l;
vars.pi_l    = pi_l;


% Diagnostics
tol = 1e-10;

% (A.66) flow consistency
res_A66 = s_l * u_l - (ml + rho) * q_l;
if abs(res_A66) > tol
    error('A.66 violated: s_l*u_l = %.3e, (ml+rho)*q_l = %.3e, residual = %.3e', ...
          s_l*u_l, (ml+rho)*q_l, res_A66);
end

% Stock identity: q_l + u_l = 1 - q_h - u_o
res_stock = q_l + u_l - (1 - q_h - u_o);
if abs(res_stock) > tol
    error('Rental stock identity violated: q_l+u_l = %.3e, RHS = %.3e, residual = %.3e', ...
          q_l+u_l, 1-q_h-u_o, res_stock);
end

% ----------------------------------------------------------
% Pre-Condition for Kappa loop: Check necessary condition given yl
% ----------------------------------------------------------

lhs_A69 = (vo * pi_h + rho) * (1 - psi) * theta_o * u_o;
rhs_A69 = m_h * q_h;

if lhs_A69 <= rhs_A69
    vars.Ry = +1e12;   % force positive residual
    return
end


% ----------------------------------------------------------
% Define n(kappa) from (A.71) and residual from (A.72)
% ----------------------------------------------------------

% helpful constant in (A.71)
numA71 = (vo*pi_h + rho) * (1-psi) * theta_o * u_o - m_h*q_h;

Z_from_kappa  = @(kap) exp(mu + sigma * norminv_noTB(kap));
Kbar_from_kap = @(kap) Kbar_from_kappa(Z_from_kappa(kap), kap, mu, sigma);

n_from_kappa = @(kap) (1/rho) * ( (numA71 ./ max(kap, realmin)) - xi*ml*q_l );

resA72 = @(kap) ...
    ( (G ./ n_from_kappa(kap)) ...
    - Fl ...
    + (1 - omega_l) * v_l * Sigma_l ) / (r + rho) ...
  + kap * ( Z_from_kappa(kap) - Kbar_from_kap(kap) ) ...
  - E;



% ----------------------------------------------------------
% Root-find κ in [0,1] using sign change of (A.72)
% ----------------------------------------------------------

epsk = 1e-10;
k0 = epsk;
k1 = 1 - epsk;

f0 = resA72(k0);
f1 = resA72(k1);

if f0 >= 0
    kappa = k0;
elseif f1 <= 0
    kappa = k1;
else
    kappa = fzero(resA72, [k0 k1]);
end

% ----------------------------------------------------------
% Back out Z, Kbar, n
% ----------------------------------------------------------
Z    = Z_from_kappa(kappa);     % (A.73)
Kbar = Kbar_from_kappa(Z, kappa, mu, sigma);   % (42) truncated mean
n    = n_from_kappa(kappa);     % (A.71)

%% ----------------------------------------------------------
% Section: Solving for the transaction threshold yl
%% ----------------------------------------------------------

% Residual yl from equation A.62
Ry = yl ...
   - Dl ...
   + Fl ...
   - (r + ml + rho)*(Cl + Cw) ...
   + xi*ml*kappa*(Z - Kbar) ...
   - (1 - omega_l + omega_l*theta_l)*v_l*Sigma_l;


%% ----------------------------------------------------------
% Return
%% ----------------------------------------------------------

vars.Ry = Ry;
vars.n = n;
vars.Kbar = Kbar;
vars.kappa = kappa;
vars.Z = Z;

end

function p = normcdf_noTB(x)
% Standard normal CDF without Statistics Toolbox
    p = 0.5 * (1 + erf(x / sqrt(2)));
end

function x = norminv_noTB(p)
% Standard normal inverse CDF without Statistics Toolbox
%
% p must be strictly between 0 and 1

    % safety clamp (important for root-finding)
    eps = 1e-12;
    p = min(max(p, eps), 1 - eps);

    x = sqrt(2) * erfinv(2 * p - 1);
end

function Kbar = Kbar_from_kappa(Z, kappa, mu, sigma)
% Average credit cost conditional on paying credit cost
% Implements equation (42) in Appendix A.3.3

% Safety for numerical stability
epsk = 1e-12;
kappa = max(kappa, epsk);

% Argument of the normal CDF
x = (log(Z) - mu - sigma^2) / sigma;

% Truncated mean of lognormal
Kbar = exp(mu + 0.5 * sigma^2) * normcdf_noTB(x) / kappa;

end
