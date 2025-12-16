function vars = a31_ownership_mkt(psi, theta_o, targets, params, verbose)
if nargin < 5, verbose = false; end
%========================================================================%
% A.3.1 Ownership-market thresholds
%========================================================================%

%% ----------------------------------------------------------
% 1. Unpack parameters
%% ----------------------------------------------------------

Fh     = params.Fh;
Fi     = params.Fi;
delta  = params.delta_h;
alpha  = params.alpha_h;
r      = params.r;
rho    = params.rho;
zeta_h = params.zeta_h;

omega_h = params.omega_h;   % seller bargaining power vs home-buyer
omega_i = params.omega_i;   % seller bargaining power vs investor

tau_h = targets.tau_h;
tau_i = targets.tau_i;

% tax-adjusted bargaining powers (Appendix A.3)
omega_h_star = omega_h / (1 + tau_h * (1 - omega_h));
omega_i_star = omega_i / (1 + tau_i * (1 - omega_i));

% meeting technology
nu_o  = params.nu_o;
eta_o = params.eta_o;
vo    = nu_o * theta_o^(-eta_o);

%% ----------------------------------------------------------
% 2. Solve for yh (equivalent to A.53)
%% ----------------------------------------------------------

y_low  = zeta_h;
y_high = zeta_h + 1e6;

yh = fzero(@yh_eq, [y_low, y_high]); % <-- must come BEFORE any printing

pi_h = (params.zeta_h / yh)^params.lambda_h;
Lambda_h_impl = 1 / pi_h;

if verbose
    fprintf('A.3.1 diagnostic:\n');
    fprintf('  yh = %.4e\n', yh);
    fprintf('  pi_h = %.6e\n', pi_h);
    fprintf('  implied Lambda_h = %.2f (target %.2f)\n', ...
            Lambda_h_impl, targets.Lambda_h);
end

%% ----------------------------------------------------------
% 3. Other objects
%% ----------------------------------------------------------

xh      = compute_xh(yh);
Sigma_h = compute_Sigma_h(yh, xh);
Sigma_i = Fi / ((1 - omega_i_star) * vo);

Ey = delta * yh;
Uh = (yh + alpha * Ey) / (r + rho + alpha);

%% ----------------------------------------------------------
% 4. Return
%% ----------------------------------------------------------

vars = struct();
vars.yh      = yh;
vars.xh      = xh;
vars.Sigma_h = Sigma_h;
vars.Sigma_i = Sigma_i;
vars.Uh      = Uh;

%% ----------------------------------------------------------
% Nested: equilibrium condition (A.50)
%% ----------------------------------------------------------
    function F = yh_eq(yh)

        xh      = compute_xh(yh);
        Sigma_h = compute_Sigma_h(yh, xh);
        Sigma_i = Fi / ((1 - omega_i_star) * vo);

        F = (xh + Fh) ...
          - (1 - omega_h_star + (1 - psi)*omega_h_star*theta_o) * vo * Sigma_h ...
          - theta_o * vo * psi * omega_i_star * Sigma_i;
    end

%% ----------------------------------------------------------
% Helper: xh(yh) from (A.52)
%% ----------------------------------------------------------
    function xh = compute_xh(yh)

        Ch = params.Ch;
        Co = params.Co;
        D  = params.D;

        den = 1 + tau_h * ((1 - psi) * omega_h_star * theta_o) / ...
                   (1 - omega_h_star + (1 - psi) * omega_h_star * theta_o) ...
                   * (r + rho + alpha)/r;

        num = yh ...
            - (r + rho + alpha) * ( ...
                Ch + (1 + tau_h)*Co - tau_h * D / r ...
                + tau_h * theta_o * omega_h_star / ...
                  (1 - omega_h_star + (1 - psi)*omega_h_star*theta_o) * ...
                  ((1 - psi)*Fh/r ...
                   + psi*(1 - omega_h_star)*omega_i_star*Fi / ...
                     (omega_h_star*(1 - omega_i_star)*r)) );

        xh = num / den;
    end

%% ----------------------------------------------------------
% Helper: Σ_h from (A.51)
%% ----------------------------------------------------------
    function Sigma_h = compute_Sigma_h(yh, xh)

        lambda_h = params.lambda_h;
        delta_h  = params.delta_h;

        coeff = zeta_h^lambda_h / ...
                ((1 + tau_h * omega_h_star) * ...
                 (r + rho + alpha) * (lambda_h - 1));

        term1 = yh^(1 - lambda_h);
        term2 = alpha * delta_h^lambda_h * xh^(1 - lambda_h) / ...
                (r + rho + alpha * (1 - delta_h^lambda_h));

        Sigma_h = coeff * (term1 + term2);
    end

end
