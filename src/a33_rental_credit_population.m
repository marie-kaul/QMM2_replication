function vars = a33_rental_credit_population(vars_A31, vars_A32, theta_o, psi, targets, params)
%========================================================================%
% Appendix A.3.3:
%   Rental-market thresholds (yl, xl),
%   Credit-cost threshold (Z),
%   City population (n),
%   Rental-market value functions Bl, Ul, Jbar,
%   Rental-market tightness theta_l (solved from investor free entry).
%
% Inputs:
%   yh, xh, Sigma_h, Sigma_i, theta_o, vo = from A.3.1–A.3.2
%   Uh      = value of owned unit (needed in Ul)
%   targets, params
%
% Outputs (in vars):
%   theta_l, yl, xl, Sigma_l, Bl, Ul, Jbar, Z, n
%========================================================================%

%% ----------------------------------------------------------
% 0. parameters and variables
% -----------------------------------------------------------

% From A.3.1 (we only need Sigma_i, Uh here; yh,xh are not used in A.3.3)
Sigma_i = vars_A31.Sigma_i;
Uh      = vars_A31.Uh; %#ok<NASGU> % may be used later if you extend to values

% From A.3.2
Sigma_h = vars_A32.Sigma_h;
pi_h    = vars_A32.pi_h;
m_h     = vars_A32.m_h;
q_h     = vars_A32.q_h;
u_o     = vars_A32.u_o;
vo      = vars_A32.vo;
G       = vars_A32.G;

% Discounting and exits
r      = params.r;
rho    = params.rho;
rho_l  = params.rho_l;

% Rental match quality distribution
zeta_l   = params.zeta_l;
lambda_l = params.lambda_l;

% Rental shocks
alpha_l = params.alpha_l;
ml      = alpha_l + rho_l;   % moving rate in rental market

% Rental costs
Dl = params.Dl;    % extra maintenance/management for landlords
Cl = params.Cl;    % landlords' transaction costs
Cw = params.Cw;    % tenant transaction costs (excl agreement fee)
Fl = params.Fl;    % flow search cost of tenants

% Ownership-side costs for (A.63)
D   = params.D;    % property maintenance cost
Co  = params.Co;   % sellers' transaction costs
Ci  = params.Ci;   % buyers' transaction costs (ownership)

% Meeting function parameters (rental)
nu_l  = params.nu_l;
eta_l = params.eta_l;

% Bargaining powers
omega_l = params.omega_l;
omega_h = params.omega_h;  % used as ω_h^* in (A.63)
omega_i = params.omega_i;  % used as ω_i^* in (A.63)

% Credit-cost distribution
xi   = params.xi;
mu   = params.mu;
sig  = params.sigma;

% Entry cost
E = params.E;

% Taxes
tau_i = targets.tau_i;      % transaction tax rate 

%% ----------------------------------------------------------
% 1. Define nested helper functions
% -----------------------------------------------------------

    function [pi_l, Sigma_l] = compute_piSigma_l(yl)
        % Step 2: π_l and Σ_l given yl (A.29)
        % π_l = (ζ_l / y_l)^{λ_l}
        pi_l = (zeta_l ./ yl).^lambda_l;
        % Σ_l = π_l y_l / [(λ_l - 1)(r + ρ + m_l)]
        Sigma_l = (pi_l .* yl) ./ ((lambda_l - 1) * (r + rho + ml));
    end

    function [s_l, theta_l, v_l] = compute_s_theta_v_l(yl, pi_l, Sigma_l)
        % Step 3: s_l, θ_l, v_l (A.63, A.64, meeting function)

        % RHS of (A.63):
        % ω_l θ_l v_l Σ_l =
        %   [1 + τ_i (1+ρ_l/r)] θ_o v_o [(1-ψ) ω_h Σ_h + ψ ω_i Σ_i]
        % + (r+ρ_l) [(1+τ_i) Co + Ci + (1+τ_i ω_i) Σ_i]
        % - τ_i (1+ρ_l/r) D

        term1 = (1 + tau_i * (1 + rho_l / r)) ...
              * theta_o * vo * ((1 - psi) * omega_h * Sigma_h + psi * omega_i * Sigma_i);

        term2 = (r + rho_l) * ( (1 + tau_i) * Co + Ci + (1 + tau_i * omega_i) * Sigma_i );

        term3 = - tau_i * (1 + rho_l / r) * D;

        RHS = term1 + term2 + term3;   % = ω_l θ_l v_l Σ_l

        % (A.64): 
        % s_l = (λ_l - 1)(r+ρ+m_l) / (ω_l y_l) * (ω_l θ_l v_l Σ_l)
        s_l = ((lambda_l - 1) * (r + rho + ml) / (omega_l * yl)) * RHS;

        % Consistency with matching: s_l = θ_l v_l π_l = ν_l π_l θ_l^{1-η_l}
        % => θ_l = (s_l / (ν_l π_l))^{1/(1-η_l)}
        theta_l = ( s_l / (nu_l * pi_l) ) .^ ( 1 / (1 - eta_l) );

        % v_l = ν_l θ_l^{-η_l}
        v_l = nu_l * theta_l .^ (-eta_l);
    end

    function [q_l, u_l] = compute_q_l_u_l(s_l)
        % Step 4: q_l, u_l from steady state (A.66)–(A.68)
        %
        % u_l = (1 - q_h - u_o) / [1 + s_l/(m_l + ρ)]
        % q_l = (s_l / (m_l + ρ)) u_l

        denom = 1 + s_l / (ml + rho);
        u_l   = (1 - q_h - u_o) / denom;
        q_l   = (s_l / (ml + rho)) * u_l;
    end

    function n_val = compute_n(kappa, q_l)
        % Step 5: city population n(κ) (A.71)
        %
        % n(κ) = (1/ρ) [ ((v_o π_h + ρ)(1-ψ) θ_o u_o - m_h q_h)/κ - ξ m_l q_l ]

        numerator = (vo * pi_h + rho) * (1 - psi) * theta_o * u_o - m_h * q_h;
        n_val = (1 / rho) * ( numerator ./ kappa - xi * ml * q_l );
    end

    function [Z, Kbar, kappa_Z_minus_Kbar] = compute_Z_Kbar(kappa)
        % Step 6: Z(κ), Kbar(κ), κ(Z-Kbar) for lognormal credit-costs
        %
        % Z = exp(mu + sig * Φ^{-1}(κ))
        % Kbar = E[K | K ≤ Z]
        % κ(Z - Kbar) = ∫_0^Z (Z - K) dΓ(K)

        eps_k = 1e-10;
        kappa_clipped = max(min(kappa, 1 - eps_k), eps_k);

        a = norminv_approx(kappa_clipped);          % Φ^{-1}(κ)
        Z = exp(mu + sig * a);               % (A.73)

        % E[K 1_{K≤Z}] for lognormal ~ LogN(mu, sig^2)
        truncated = exp(mu + 0.5 * sig^2) * normcdf_approx(a - sig);

        % Conditional mean
        Kbar = truncated ./ kappa_clipped;

        kappa_Z_minus_Kbar = kappa .* (Z - Kbar);
    end

    function [kappa_star, Z_star, Kbar_star, n_star] = solve_kappa(yl, Sigma_l, v_l, q_l)
        % Step 7: solve city-indifference condition (A.72) for κ
        %
        % (G/n) - F_l + (1-ω_l)v_l Σ_l/(r+ρ) + κ(Z-Kbar) - E = 0

        rental_term = (1 - omega_l) * v_l * Sigma_l / (r + rho);

        % residual F_kappa(κ)
        function Fk = kappa_residual(kappa)
            n_val = compute_n(kappa, q_l);
            [Z, Kbar, kZminusKbar] = compute_Z_Kbar(kappa);
            % Note kZminusKbar = κ(Z - Kbar)
            Fk = (G / n_val) - Fl + rental_term + kZminusKbar - E;
        end

        % bracket κ in (0,1)
        k_low  = 1e-6;
        k_high = 1 - 1e-6;

        F_low  = kappa_residual(k_low);
        F_high = kappa_residual(k_high);

       if sign(F_low) == sign(F_high)
            % κ is infeasible for this yl
            % Return a flag to yl-equilibrium that this yl is invalid
            kappa_star = NaN;
            Z_star = NaN;
            Kbar_star = NaN;
            n_star = NaN;
            return;
        end


        kappa_star = fzero(@kappa_residual, [k_low, k_high]);

        % Once κ* is found, compute Z*, Kbar*, n*
        n_star = compute_n(kappa_star, q_l);
        [Z_star, Kbar_star, ~] = compute_Z_Kbar(kappa_star);
    end

    function [Fyl, pi_l, Sigma_l, s_l, theta_l, v_l, q_l, u_l, ...
              kappa, Z, Kbar, n_val] = yl_equilibrium(yl)
        % Full equilibrium map for a given yl:
        %  - computes all rental and credit vars
        %  - returns F(yl) (A.62)

        % Step 2: π_l and Σ_l
        [pi_l, Sigma_l] = compute_piSigma_l(yl);

        % Step 3: s_l, θ_l, v_l
        [s_l, theta_l, v_l] = compute_s_theta_v_l(yl, pi_l, Sigma_l);

        % Step 4: q_l, u_l
        [q_l, u_l] = compute_q_l_u_l(s_l);

        % Step 7: κ, Z, Kbar, n
        [kappa, Z, Kbar, n_val] = solve_kappa(yl, Sigma_l, v_l, q_l);

        if isnan(kappa)
            % yl is infeasible; return a very large positive Fyl
            Fyl = +1e9 + (yl - zeta_l)*1e3;
            return;
        end

        % Step 8: threshold residual F(yl) (A.62)
        %
        % F(yl) = yl - D_l + F_l
        %         - (r + m_l + ρ)(C_l + C_w)
        %         + ξ m_l κ(Z - Kbar)
        %         - (1 - ω_l + ω_l θ_l) v_l Σ_l

        direct_term = yl - Dl + Fl ...
                      - (r + ml + rho) * (Cl + Cw);

        credit_term = xi * ml * kappa * (Z - Kbar);

        surplus_term = (1 - omega_l + omega_l * theta_l) * v_l * Sigma_l;

        Fyl = direct_term + credit_term - surplus_term;
    end

    function F = yl_residual(yl)
        F = yl_equilibrium(yl);
    end


% %% ----------------------------------------------------------
% % 2. Solve for yl using fzero on F(yl)
% % -----------------------------------------------------------
% 
% % Start extremely close to zeta_l
% yl_lower = zeta_l * 1.0001;
% F_low = yl_residual(yl_lower);
% 
% % Initialize upper bound
% yl_upper = yl_lower * 2;
% F_high = yl_residual(yl_upper);
% 
% % Expand until sign changes or max level hit
% max_expand = 60;   % MUCH larger than before!
% iter = 0;
% 
% while sign(F_low) == sign(F_high) && iter < max_expand
%     yl_upper = yl_upper * 2;
%     F_high = yl_residual(yl_upper);
%     iter = iter + 1;
% end
% 
% if iter == max_expand
%     fprintf('\n\n*** DEBUG: yl never bracketed within expansion ***\n');
%     fprintf('yl_lower = %.4e, F_low = %.4e\n', yl_lower, F_low);
%     fprintf('yl_upper = %.4e, F_high = %.4e\n', yl_upper, F_high);
%     error('A.3.3: Could NOT bracket yl. Increase max_expand or inspect Fyl.');
% end
% 
% fprintf('yl_lower = %.4e, F_low = %.4e\n', yl_lower, F_low);
% fprintf('yl_upper = %.4e, F_high = %.4e\n', yl_upper, F_high);
% 
% 
% % Now solve with fzero
% yl_star = fzero(@yl_residual, [yl_lower, yl_upper]);

% ----------------------------------------------------------
% Solve for yl using bisection on [yl_lower, yl_upper]
% ----------------------------------------------------------

% Start extremely close to zeta_l
yl_lower = zeta_l * 1.0001;
F_low = yl_residual(yl_lower);

% Initialize upper bound
yl_upper = yl_lower * 2;
F_high = yl_residual(yl_upper);

max_iter = 100;
tol_F    = 1e-6;
tol_y    = 1e-6;

yl_low  = yl_lower;
yl_high = yl_upper;
F_low   = yl_residual(yl_low);
F_high  = yl_residual(yl_high);

if F_low * F_high > 0
    error('A.3.3: Bisection called without a valid sign change.');
end

for it = 1:max_iter
    yl_mid = 0.5 * (yl_low + yl_high);
    F_mid  = yl_residual(yl_mid);

    % Check convergence
    if abs(F_mid) < tol_F || abs(yl_high - yl_low) < tol_y
        yl_star = yl_mid;
        break;
    end

    % Narrow the bracket
    if F_low * F_mid < 0
        yl_high = yl_mid;
        F_high = F_mid;
    else
        yl_low = yl_mid;
        F_low  = F_mid;
    end

    if it == max_iter
        warning('A.3.3: Bisection did not fully converge, using last yl_mid.');
        yl_star = yl_mid;
    end
end


%% ----------------------------------------------------------
% 3. Recompute all equilibrium objects at yl*
% -----------------------------------------------------------

[F_yl, pi_l_star, Sigma_l_star, s_l_star, theta_l_star, v_l_star, ...
 q_l_star, u_l_star, kappa_star, Z_star, Kbar_star, n_star] = yl_equilibrium(yl_star); %#ok<NASGU>


%% ----------------------------------------------------------
% 4. Pack results
% -----------------------------------------------------------

vars = struct();

vars.yl       = yl_star;
vars.kappa    = kappa_star;
vars.Z        = Z_star;
vars.Kbar     = Kbar_star;
vars.n        = n_star;

vars.theta_l  = theta_l_star;
vars.Sigma_l  = Sigma_l_star;
vars.v_l      = v_l_star;
vars.s_l      = s_l_star;
vars.q_l      = q_l_star;
vars.u_l      = u_l_star;

end

%% ----------------------------------------------------------
% Helpers
% -----------------------------------------------------------

function x = norminv_approx(p)
    % Approximation of the inverse normal CDF 
    % Accuracy ~ 1e-12. Valid for 0 < p < 1.

    % Coefficients from Wichura (Algorithm AS 241)
    a = [ -3.969683028665376e+01,  2.209460984245205e+02,...
          -2.759285104469687e+02,  1.383577518672690e+02,...
          -3.066479806614716e+01,  2.506628277459239e+00 ];

    b = [ -5.447609879822406e+01,  1.615858368580409e+02,...
          -1.556989798598866e+02,  6.680131188771972e+01,...
          -1.328068155288572e+01 ];

    c = [ -7.784894002430293e-03, -3.223964580411365e-01,...
          -2.400758277161838e+00, -2.549732539343734e+00,...
           4.374664141464968e+00,  2.938163982698783e+00 ];

    d = [ 7.784695709041462e-03,  3.224671290700398e-01,...
          2.445134137142996e+00,  3.754408661907416e+00 ];

    % Define break-points
    plow = 0.02425;
    phigh = 1 - plow;

    if p < plow
        q = sqrt(-2*log(p));
        x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / ...
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1);
    elseif p > phigh
        q = sqrt(-2*log(1-p));
        x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / ...
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1);
    else
        q = p - 0.5;
        r = q*q;
        x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / ...
            (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1);
    end
end

function y = normcdf_approx(x)
    % Fast approximation of normal CDF (Abramowitz-Stegun)
    % Accuracy ~ 1e-7 (more than enough for calibration).

    t = 1 ./ (1 + 0.2316419 * abs(x));
    d = 0.3989423 * exp(-x .* x / 2);   % standard normal PDF

    poly = (((1.330274429 * t ...
            - 1.821255978) .* t ...
            + 1.781477937) .* t ...
            - 0.356563782) .* t ...
            + 0.319381530;

    y = 1 - d .* poly;

    y(x < 0) = 1 - y(x < 0);
end
