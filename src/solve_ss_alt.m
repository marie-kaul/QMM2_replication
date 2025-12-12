%========================================================================%
% This file solves for steady state.
% Input: data/calibration_targets.mat.
% Output: 
% 
% From Section A.3 of the paper: 
% The solution method is based on a numerical search over the fraction ψ of 
% investors among buyers and ownership-market tightness θo that satisfy two 
% equations representing equilibrium in the ownership and rental markets. 
% Within this search, given a (ψ,θo), the ownership-market thresholds (yh, xh) 
% and rental-market and credit-cost thresholds (yl ,Z) are found by solving 
% two equations numerically.
%========================================================================%


clear; clc;

load('../data/calibration_data.mat');

%% ----------------------------------------------------------
% Initial guesses (interior!)
%% ----------------------------------------------------------

psi0     = 0.1;     % much larger investor share
theta_o0 = 3.0;     % tighter ownership market

% transform to unconstrained variables
y0 = zeros(2,1);
y0(1) = log(psi0/(1-psi0));   % logit transform
y0(2) = log(theta_o0);        % log transform

%% ----------------------------------------------------------
% Solver options
%% ----------------------------------------------------------

opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',3000);

%% ----------------------------------------------------------
% Solve for (psi, theta_o)
%% ----------------------------------------------------------

[y_star, ~, exitflag, output] = fsolve(@outer_residual, y0, opts);

if exitflag <= 0
    error('Equilibrium solver did not converge');
end

%% ----------------------------------------------------------
% Recover equilibrium values
%% ----------------------------------------------------------

psi     = 1/(1 + exp(-y_star(1)));
theta_o = exp(y_star(2));

fprintf('\n==== STEADY STATE FOUND ====\n');
fprintf('psi      = %.6f\n', psi);
fprintf('theta_o = %.6f\n', theta_o);

%% ----------------------------------------------------------
% Compute final equilibrium objects
%% ----------------------------------------------------------

[vars_A31, vars_A32, vars_A33, resA34] = evaluate_all(psi, theta_o);

fprintf('\nResiduals:\n');
fprintf('A.35 = %.3e\n', resA34(1));
fprintf('A.74 = %.3e\n', resA34(2));

%% ----------------------------------------------------------
% Save equilibrium
%% ----------------------------------------------------------

steady = struct();
steady.psi      = psi;
steady.theta_o  = theta_o;
steady.vars_A31 = vars_A31;
steady.vars_A32 = vars_A32;
steady.vars_A33 = vars_A33;
steady.resA34   = resA34;

save('../data/steady_state.mat','steady');

disp('Steady state saved to data/steady_state.mat');

%% ==========================================================
% Nested functions
%% ==========================================================

function F = outer_residual(y)

    % map unconstrained -> economic variables
    psi_i     = 1/(1 + exp(-y(1)));
    theta_o_i = exp(y(2));

    % safety clamps
    psi_i     = min(max(psi_i,1e-10),1-1e-10);
    theta_o_i = max(theta_o_i,1e-10);

    try
        [~, ~, ~, resA34_i] = evaluate_all(psi_i, theta_o_i);
    catch
        % if A.3.3 fails (e.g. no unique yl), penalize
        F = [1e6; 1e6];
        return
    end

    if any(~isfinite(resA34_i))
        F = [1e6; 1e6];
    else
        F = resA34_i;
    end
end


function [vars_A31, vars_A32, vars_A33, resA34] = evaluate_all(psi, theta_o)

    % A.3.1 Ownership thresholds
    vars_A31 = a31_ownership_mkt(psi, theta_o, targets, params);

    % A.3.2 Ownership stocks and flows
    vars_A32 = a32_ownership_other(vars_A31, psi, theta_o, targets, params);

    % A.3.3 Rental market, credit, population
    vars_A33 = a33_rental_credit_population1( ...
        vars_A31, vars_A32, theta_o, psi, targets, params);

    % A.3.4 residuals
    resA34 = a34_eq_target(vars_A33, vars_A32, theta_o, psi, targets, params);

end


















