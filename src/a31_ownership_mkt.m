function [vars] = a31_ownership_mkt(psi, theta_o, targets)

%========================================================================%
% A.3.1 Ownership-market thresholds

% Input: (ψ,θo), targets (exogenously calibrated vars)
% Outputs: ...

% This part of the solution method derives an equation satisfied by the ownership-market transaction
% threshold yh, which can be solved taking as given (ψ,θo). Once yh is known, the moving threshold
% xh is determined, along with other variables related to the ownership market
%========================================================================%

%% 1. Unpack needed parameters from targets
r    = targets.rg;        % discount rate r
rho  = targets.rho;       % depopulation rate ρ
alpha_h = targets.alpha_h;
lambda_h = targets.lambda_h;
delta_h  = targets.delta_h;
Fh   = targets.Fh;
Fi   = targets.Fi;
Ch   = targets.Ch;
Co   = targets.Co;
D    = targets.D_over_P;    % or D, depending on your units
tau_h = targets.tau_h;
zeta_h = targets.zeta_h;    % lower bound for yh
vo   = targets.vo;          % meeting rate vo(θo) = υo θo^{-ηo}
uo   = targets.uo;          % stock of sellers, needed later
omega_h = targets.omega_star_h;
omega_star_i = targets.omega_star_i;
nu_o  = targets.uo; % (if needed for stock-flow)

%% 2. Solve for y_h
y_low = zeta_h;
y_high = 50;   % pick something large: LHS→positive for large yh

yh = fzero(@yh_eq, [y_low, y_high]);    % function yh_eq defined below

%% 4. compute xh, pi_h, Sigma_h, Sigma_i as outputs
xh = compute_xh(yh, psi, theta_o, targets);
pi_h = (zeta_h / yh)^lambda_h;
Sigma_h = compute_Sigma_h(yh, xh, targets);
Sigma_i = Fi / ((1 - omega_star_i)*vo);

%% 5. Pack results 
vars = struct();
vars.yh = yh;
vars.xh = xh;
vars.pi_h = pi_h;
vars.Sigma_h = Sigma_h;
vars.Sigma_i = Sigma_i;

%% ----------------------------------------------------------
% Helper functions
%% ----------------------------------------------------------

% the root-finding equation in yh (implementing A.53)

    function F = yh_eq(yh)

        % ---- Compute xh from (A.52) ---- %
        xh = compute_xh(yh, psi, theta_o, targets);     % compute xh defined below

        % ---- Compute Σ_h(yh,xh) from (A.51) ---- %
        Sigma_h = compute_Sigma_h(yh, xh, targets);     % compute sigma h defined below

        % ---- Σ_i from (38): Σi = Fi /((1 - ω_i*) vo) ---- %
        Sigma_i = Fi / ((1 - omega_star_i)*vo);

        % ---- Left-hand side of (A.53) ---- %
        LHS = xh + Fh ...
              - (1 - omega_h + (1-psi)*omega_h*theta_o)*vo * Sigma_h ...
              - psi*theta_o * omega_star_i/(1 - omega_star_i) * Fi;

        F = LHS;
    end

% compute xh from A.52

function xh = compute_xh(yh, psi, theta_o, T)
    r = T.rg; rho = T.rho; alpha_h = T.alpha_h;
    Ch = T.Ch; Co = T.Co; tau_h = T.tau_h;
    D = T.D_over_P; Fh = T.Fh; Fi = T.Fi;
    omega_h = T.omega_star_h; omega_star_i = T.omega_star_i;
    vo = T.vo;

    num = yh - (r+rho+alpha_h)*( ...
            Ch + (1+tau_h)*Co - tau_h*D/r ...
            + tau_h * theta_o*omega_h /(1 - omega_h + (1-psi)*omega_h*theta_o) * ...
              ((1-psi)*Fh/r + psi*(1-omega_h)*omega_star_i*Fi/(omega_h*(1-omega_star_i)*r)));

    den = 1 + tau_h*((1-psi)*omega_h*theta_o)/(1 - omega_h + (1-psi)*omega_h*theta_o) ...
            * (r+rho+alpha_h)/r;

    xh = num/den;
end

% compute sigma h from A 51

function Sigma_h = compute_Sigma_h(yh, xh, T)
    r=T.rg; rho=T.rho; alpha_h=T.alpha_h;
    lambda_h=T.lambda_h; delta_h=T.delta_h; zeta_h=T.zeta_h;
    omega_h=T.omega_star_h; tau_h=T.tau_h; vo=T.vo; theta_o=T.theta_o;

    coeff = zeta_h^lambda_h / ((1 + tau_h*omega_h)*(r+rho+alpha_h)*(lambda_h-1));
    term1 = yh^(1 - lambda_h);
    term2 = alpha_h * delta_h^lambda_h * xh^(1 - lambda_h) / (r+rho + alpha_h*(1 - delta_h^lambda_h));
    Sigma_h = coeff * vo * (term1 + term2);
end


end