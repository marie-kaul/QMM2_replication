function vars = a31_ownership_mkt(psi, theta_o, targets, params)
%========================================================================%
% A.3.1 Ownership-market thresholds
%
% Inputs:
%   psi,theta_o = guesses for equilibrium vars
%   targets     = Table 3 calibration targets 
%   params      = Table 4 calibrated parameters 
%
% Outputs:
% vars.yh = yh;               % ownership-market transaction threshold
% vars.xh = xh;               % moving threshold
% vars.Sigma_i = Sigma_i;     % seller share of surplus when selling to investor
%
%========================================================================%

%% ----------------------------------------------------------
% 1. Unpack calibrated parameters 
%% ----------------------------------------------------------

Fh      = params.Fh;
Fi      = params.Fi;
delta   = params.delta_h;
alpha   = params.alpha_h;
r       = params.r;
rho     = params.rho;
zeta_h  = params.zeta_h;

omega_h     = params.omega_h;       % bargaining power (home-buyer)
omega_star_i = params.omega_i;      % bargaining power (investor)

% Meeting function parameters: vo = ν_o * θ_o^{−η_o}
nu_o   = params.nu_o;
eta_o  = params.eta_o;
vo     = nu_o * theta_o^(-eta_o);

%% ----------------------------------------------------------
% 2. Solve for yh using (A.53)
%% ----------------------------------------------------------

y_low  = zeta_h;
y_high = zeta_h + 1e6;   % large enough upper bound

yh = fzero(@yh_eq, [y_low, y_high]); % solver for nonlinear eq

%% ----------------------------------------------------------
% 3. Compute xh, pi_h, Sigma_h, Sigma_i, Uh
%% ----------------------------------------------------------

xh       = compute_xh(yh, psi, theta_o, targets, params, vo);
Sigma_h  = compute_Sigma_h(yh, xh, targets, params, vo, theta_o);
Sigma_i  = Fi / ((1 - omega_star_i) * vo);

Ey = delta * yh;  
Uh = ( yh + alpha * Ey ) / (r + rho + alpha);

%% ----------------------------------------------------------
% 4. Return output
%% ----------------------------------------------------------

vars = struct();
vars.yh = yh;               % ownership-market transaction threshold
vars.xh = xh;               % moving threshold
vars.Sigma_i = Sigma_i;     % seller share of surplus when selling to investor
vars.Uh      = Uh;          % the value of being a homeowner

%% ----------------------------------------------------------
% Nested: equation (A.53) for root-finding of yh
%% ----------------------------------------------------------
    function F = yh_eq(yh)

        xh = compute_xh(yh, psi, theta_o, targets, params, vo);
        Sigma_h = compute_Sigma_h(yh, xh, targets, params, vo, theta_o);
        Sigma_i = Fi / ((1 - omega_star_i) * vo);

        % (A.53) LHS
        LHS = xh + Fh ...
            - (1 - omega_h + (1 - psi) * omega_h * theta_o) * vo * Sigma_h ...
            - psi * theta_o * omega_star_i/(1 - omega_star_i) * Fi;

        F = LHS;
    end

end % end main function

%==========================================================================
% Helper: compute xh from (A.52)
%==========================================================================

function xh = compute_xh(yh, psi, theta_o, T, P, vo)

    r = P.r; rho = P.rho; alpha_h = P.alpha_h;

    Ch = P.Ch; Co = P.Co; D = P.D; Fh = P.Fh; Fi = P.Fi;
    tau_h = T.tau_h;

    omega_h      = P.omega_h;
    omega_star_i = P.omega_i;

    % Denominator term: 1 + τ_h*((1−ψ)ω_h θ_o)/(1 − ω_h + (1−ψ)ω_h θ_o)*(r+ρ+α_h)/r
    den = 1 + tau_h * ((1 - psi) * omega_h * theta_o) / ...
               (1 - omega_h + (1 - psi) * omega_h * theta_o) * (r + rho + alpha_h)/r;

    % Long numerator, following (A.52)
    num = yh ...
        - (r + rho + alpha_h) * ( Ch + (1 + tau_h)*Co - tau_h * D / r ...
        + tau_h * theta_o * omega_h / (1 - omega_h + (1 - psi)*omega_h*theta_o) * ...
          ((1 - psi)*Fh/r + psi*(1 - omega_h)*omega_star_i * Fi / ...
          (omega_h * (1 - omega_star_i) * r)) );

    xh = num / den;

end

%==========================================================================
% Helper: compute Σ_h from (A.51)
%==========================================================================

function Sigma_h = compute_Sigma_h(yh, xh, targets, params, vo, theta_o)

    r = params.r; rho = params.rho; alpha_h = params.alpha_h;
    lambda_h = params.lambda_h; delta_h = params.delta_h;

    zeta_h  = params.zeta_h;
    tau_h   = targets.tau_h;
    omega_h = params.omega_h;

    coeff = zeta_h^lambda_h / ((1 + tau_h * omega_h) * (r + rho + alpha_h) * (lambda_h - 1));

    term1 = yh^(1 - lambda_h);
    term2 = alpha_h * delta_h^lambda_h * xh^(1 - lambda_h) / ...
        (r + rho + alpha_h * (1 - delta_h^lambda_h));

    Sigma_h = coeff * vo * (term1 + term2);

end
