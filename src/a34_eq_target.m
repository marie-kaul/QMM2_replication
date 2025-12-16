function res = a34_eq_target(vars_A33, vars_A32, theta_o, psi, targets, params)

%========================================================================%
% A.3.4 Solving for psi and theta_o
% Residuals:
%   (1) (A.35) relationship between market tightnesses across markets
%   (2) (A.74) marginal home-buyer indifference condition
%========================================================================%

%% Parameters
r       = params.r;
rho     = params.rho;
Fh      = params.Fh;
Fl      = params.Fl;

omega_h = params.omega_h;
omega_l = params.omega_l;

% tax-adjusted omega_h* (Appendix A.3 / A.74)
tau_h = targets.tau_h;
omega_h_star = omega_h / (1 + tau_h * (1 - omega_h));

%% Variables from previous blocks
theta_l = vars_A33.theta_l;

vo      = vars_A32.vo;
vl      = vars_A33.v_l;

Sigma_h = vars_A32.Sigma_h;
Sigma_l = vars_A33.Sigma_l;

Z       = vars_A33.Z;

u_o     = vars_A32.u_o;
u_l     = vars_A33.u_l;

n       = vars_A33.n;

%% (A.35): tightness relationship
res_A35 = (((1 - psi) * theta_o - 1) * u_o) ...
        + ((theta_l - 1) * u_l) ...
        - (n - 1);

%% (A.74): marginal buyer indifference
res_A74 = ( (1 - omega_h_star) * vo * Sigma_h ...
          - (1 - omega_l)      * vl * Sigma_l ) ...
        - ( (r + rho) * Z + Fh - Fl );

%% Collect residuals
res = [res_A35; res_A74];

end
