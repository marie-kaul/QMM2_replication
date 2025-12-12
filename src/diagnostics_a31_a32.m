function diag = diagnostics_a31_a32(vars_A31, vars_A32, psi, theta_o, targets, params)
%==========================================================================
% Diagnostics for Appendix A.3.1–A.3.2 (Ownership Market)
%
% This function prints and returns key residuals and implied objects used
% to verify correctness of the ownership-market block.
%
% Inputs:
%   vars_A31   struct from a31_ownership_mkt
%   vars_A32   struct from a32_ownership_other
%   psi        investor share
%   theta_o    ownership-market tightness
%   targets    calibration targets
%   params     structural parameters
%
% Output:
%   diag       struct collecting diagnostic values
%==========================================================================

fprintf('\n================ A.3 OWNERSHIP DIAGNOSTICS ================\n');

%% ----------------------------------------------------------
% A.3.1: Transaction probability and Lambda_h
%% ----------------------------------------------------------

pi_h = vars_A32.pi_h;
Lambda_h_impl = 1 / pi_h;
Lambda_gap = abs(Lambda_h_impl - targets.Lambda_h);

fprintf('pi_h = %.6f\n', pi_h);
fprintf('Lambda_h = %.4f (target %.4f)\n', ...
        Lambda_h_impl, targets.Lambda_h);
fprintf('Gap in Lambda_h = %.4f\n', Lambda_gap);

if Lambda_gap > 0.5
    warning('Lambda_h gap is large (%.4f). Expected away from equilibrium.', Lambda_gap);
end

%% ----------------------------------------------------------
% A.3.1: Threshold ordering
%% ----------------------------------------------------------

delta_yh = params.delta_h * vars_A31.yh;
fprintf('Check δ_h y_h < x_h : %.4f < %.4f\n', ...
        delta_yh, vars_A31.xh);

%% ----------------------------------------------------------
% A.3.1: Equation (A.50) residual
%% ----------------------------------------------------------

lhs_A50 = vars_A31.xh + params.Fh;

rhs_A50 = (1 - params.omega_h + (1 - psi)*params.omega_h*theta_o) ...
          * vars_A32.vo * vars_A32.Sigma_h ...
          + theta_o * vars_A32.vo * psi * params.omega_i * vars_A31.Sigma_i;

A50_resid = lhs_A50 - rhs_A50;
fprintf('A.50 residual (should be ~0): %.3e\n', A50_resid);

%% ----------------------------------------------------------
% A.3.2: Moving hazard m_h (A.58)
%% ----------------------------------------------------------

fprintf('m_h = %.6f (target 1/Tmh = %.6f)\n', ...
        vars_A32.m_h, 1 / targets.Tmh);

%% ----------------------------------------------------------
% A.3.2: Sales rate consistency (eq. 36)
%% ----------------------------------------------------------

s_o_check = theta_o * vars_A32.vo * (psi + (1 - psi)*pi_h);

fprintf('s_o consistency | model %.6f, implied %.6f\n', ...
        vars_A32.s_o, s_o_check);

%% ----------------------------------------------------------
% A.3.2: Stock-flow consistency (A.56)
%% ----------------------------------------------------------

lhs_A56 = (1 - psi) * vars_A32.s_o * vars_A32.u_o;
rhs_A56 = (vars_A32.m_h + params.rho) * vars_A32.q_h;

A56_resid = lhs_A56 - rhs_A56;
fprintf('A.56 residual (should be ~0): %.3e\n', A56_resid);

%% ----------------------------------------------------------
% A.3.2: Mass sanity
%% ----------------------------------------------------------

fprintf('u_o + q_h = %.6f (<= 1)\n', ...
        vars_A32.u_o + vars_A32.q_h);

%% ----------------------------------------------------------
% Prices sanity (A.54)–(A.55)
%% ----------------------------------------------------------

r  = params.r;
Co = params.Co;
D  = params.D;

omega_h = params.omega_h;
omega_i = params.omega_i;

Sigma_h = vars_A32.Sigma_h;
Sigma_i = vars_A31.Sigma_i;
vo      = vars_A32.vo;

Ph = ((r + theta_o*vo*(1-psi)*pi_h)/r) * (omega_h*Sigma_h/pi_h) ...
     + (theta_o*vo*psi*omega_i*Sigma_i)/r ...
     + Co - D/r;

Pi = Co ...
     + (theta_o*vo*((1-psi)*omega_h*Sigma_h + psi*omega_i*Sigma_i) - D)/r ...
     + omega_i*Sigma_i;

fprintf('Ph = %.6e (should be >0)\n', Ph);
fprintf('Pi = %.6e (should be >0)\n', Pi);

%% ----------------------------------------------------------
% Sales decomposition
%% ----------------------------------------------------------

S_total = vars_A32.s_o * vars_A32.u_o;
i_share = psi / (psi + (1-psi)*pi_h);

fprintf('S_total = %.6e\n', S_total);
fprintf('Investor share i = %.6e (in (0,1))\n', i_share);

fprintf('===========================================================\n');

%% ----------------------------------------------------------
% Collect diagnostics
%% ----------------------------------------------------------

diag = struct();
diag.pi_h          = pi_h;
diag.Lambda_h      = Lambda_h_impl;
diag.Lambda_gap    = Lambda_gap;
diag.A50_resid     = A50_resid;
diag.A56_resid     = A56_resid;
diag.m_h           = vars_A32.m_h;
diag.s_o_resid     = vars_A32.s_o - s_o_check;
diag.Ph            = Ph;
diag.Pi            = Pi;
diag.S_total       = S_total;
diag.i_share       = i_share;

end
