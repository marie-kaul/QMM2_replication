%========================================================================%
% This file solves for steady state (outer loop over psi, theta_o)
% using fminsearch on the A.3.4 equilibrium conditions.
%========================================================================%

clear; clc;

load('/Users/ewaste/Documents/PhD/Replication audit/data/calibration_data.mat');

%% ----------------------------------------------------------
% Initial guesses (interior)
%% ----------------------------------------------------------
psi0     = 0.003;
theta0   = 1.0;

% Unconstrained parameterization
y0 = [
    log(psi0/(1-psi0));   % logit(psi)
    log(theta0)           % log(theta_o)
];

%% ----------------------------------------------------------
% Test objective at initial guess
%% ----------------------------------------------------------
test_J = outer_objective(y0, targets, params);
fprintf('Initial objective J = %.4e\n', test_J);

%% ----------------------------------------------------------
% Outer solve with fminsearch
%% ----------------------------------------------------------
opts = optimset( ...
    'Display','iter', ...
    'MaxIter', 800, ...
    'MaxFunEvals', 4000, ...
    'TolX', 1e-10, ...
    'TolFun', 1e-10 );

[y_sol, Jstar, exitflag, output] = ...
    fminsearch(@(y) outer_objective(y, targets, params), y0, opts);

%% ----------------------------------------------------------
% Map back to economic parameters
%% ----------------------------------------------------------
psi     = 1/(1 + exp(-y_sol(1)));
theta_o = exp(y_sol(2));

% IMPORTANT: clamp again (must match objective!)
psi = min(max(psi, 1e-5), 1 - 1e-8);

fprintf('\n=== OUTER SOLUTION (fminsearch) ===\n');
fprintf('psi     = %.6e\n', psi);
fprintf('theta_o = %.8f\n', theta_o);
fprintf('J       = %.6e\n', Jstar);
fprintf('exitflag= %d\n', exitflag);

%% ----------------------------------------------------------
% Final full evaluation at solution (verbose)
%% ----------------------------------------------------------
vars_A31 = a31_ownership_mkt(psi, theta_o, targets, params, true);
vars_A32 = a32_ownership_other(vars_A31, psi, theta_o, targets, params);
vars_A33 = a33_rental_credit_population1(vars_A31, vars_A32, theta_o, psi, targets, params);
res      = a34_eq_target(vars_A33, vars_A32, theta_o, psi, targets, params);

fprintf('\n=== CHECK AT SOLUTION (unscaled) ===\n');
fprintf('A.35 residual = %.6e\n', res(1));
fprintf('A.74 residual = %.6e\n', res(2));
fprintf('n (implied)   = %.6f (target 1)\n', vars_A33.n);

%% ----------------------------------------------------------
% Inequality checks (Appendix A.3)
%% ----------------------------------------------------------
yh = vars_A31.yh;
xh = vars_A31.xh;

if params.delta_h * yh < xh
    disp('A.3.1: δ_h y_h < x_h is satisfied.')
else
    error('A.3.1: δ_h y_h < x_h is NOT satisfied.');
end

diag_A3 = diagnostics_a31_a32(vars_A31, vars_A32, psi, theta_o, targets, params);

%% ----------------------------------------------------------
% Save steady state
%% ----------------------------------------------------------
SS = struct( ...
    'psi',psi, ...
    'theta_o',theta_o, ...
    'vars_A31',vars_A31, ...
    'vars_A32',vars_A32, ...
    'vars_A33',vars_A33, ...
    'residuals',res, ...
    'Jstar',Jstar, ...
    'exitflag',exitflag, ...
    'output',output, ...
    'diagnostics',diag_A3 );

save('steady_state_solution.mat','SS');

fprintf('\nSteady state saved to steady_state_solution.mat\n');

%% ==========================================================
% SANITY CHECK: key parameters and units
% ==========================================================

fprintf('\n================ UNIT SANITY CHECK =================\n');

fprintf('--- Monetary units (should be THOUSANDS of $) ---\n');
fprintf('P_avg      = %.3f   (paper: 402)\n', targets.P_avg);
fprintf('Fh         = %.3f   (paper: 9.84)\n', params.Fh);
fprintf('Fl         = %.3f   (paper: 9.72)\n', params.Fl);
fprintf('Co         = %.3f   (paper: 18.1)\n', params.Co);
fprintf('Cl         = %.3f   (paper: 2.28)\n', params.Cl);
fprintf('Cw         = %.3f   (paper: 0.709)\n', params.Cw);
fprintf('D          = %.3f   (paper: 10.5)\n', params.D);
fprintf('Dl         = %.3f   (paper: 2.20)\n', params.Dl);
fprintf('E          = %.3f   (paper: 21.1)\n', params.E);

fprintf('\n--- Match quality thresholds ---\n');
fprintf('zeta_h     = %.3f   (paper: 33.6)\n', params.zeta_h);
fprintf('zeta_l     = %.3f   (paper: 24.6)\n', params.zeta_l);

fprintf('\n--- Credit-cost distribution (log of THOUSANDS) ---\n');
fprintf('mu         = %.3f   (paper: 5.05)\n', params.mu);
fprintf('sigma      = %.3f   (paper: 0.674)\n', params.sigma);

fprintf('\n--- Interest / exit rates ---\n');
fprintf('r          = %.4f\n', params.r);
fprintf('rho        = %.4f\n', params.rho);
fprintf('rho_l      = %.4f\n', params.rho_l);

fprintf('===================================================\n\n');


%% =======================================================================
% Objective function for fminsearch
% =======================================================================

function J = outer_objective(y, targets, params)

    % Unconstrained -> constrained
    psi     = 1/(1 + exp(-y(1)));
    theta_o = exp(y(2));

    % Hard interior bounds for numerical stability
    psi = min(max(psi, 1e-5), 1 - 1e-8);

    % Soft bounds on theta_o
    pen = 0;
    if theta_o < 1e-3
        pen = pen + 1e6 * (1e-3 - theta_o)^2;
    end
    if theta_o > 50
        pen = pen + 1e6 * (theta_o - 50)^2;
    end

    try
        vars_A31 = a31_ownership_mkt(psi, theta_o, targets, params, false);
        vars_A32 = a32_ownership_other(vars_A31, psi, theta_o, targets, params);
        vars_A33 = a33_rental_credit_population1(vars_A31, vars_A32, theta_o, psi, targets, params);

        res = a34_eq_target(vars_A33, vars_A32, theta_o, psi, targets, params);

        % Scale A.74 so it doesn't dominate
        res_scaled = [res(1); res(2)/1e4];

        if any(~isfinite(res_scaled))
            J = 1e12;
        else
            J = sum(res_scaled.^2) + pen;
        end

    catch
        J = 1e12 + pen;
    end
end
