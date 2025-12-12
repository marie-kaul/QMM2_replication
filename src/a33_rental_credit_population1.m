function vars = a33_rental_credit_population1(vars_A31, vars_A32, ...
    theta_o, psi, targets, params)
%========================================================================%
% Appendix A.3.3:
% Rental-market thresholds yl, credit-cost threshold Z,
% city population n, and rental tightness theta_l.
%
% Conceptually, we have
% 1. Outer loop over values of yl
% 2. Given yl, find kappa
% 3. Given kappa, we can determine Z and n
%========================================================================%

zeta_l   = params.zeta_l;

%% ----------------------------------------------------------
% Initialize loop over yl
%% ----------------------------------------------------------

% Initial bracket (must satisfy sign change)
yl_lo = zeta_l * 1.001;   
yl_hi = 5 * yl_lo;    

vars_lo = yl_residual(yl_lo, vars_A31, vars_A32, ...
    theta_o, psi, targets, params);
R_lo = vars_lo.Ry;

% ----------------------------------------------------------
% Uniqueness / existence check for y_l (Appendix A.3.3)
% ----------------------------------------------------------
if R_lo >= 0
    error(['A.3.3 uniqueness/existence condition failed: ', ...
           'Residual at y_l = zeta_l is non-negative (R = %.3e). ', ...
           'No unique y_l exists.'], R_lo);
end


vars_hi = yl_residual(yl_hi, vars_A31, vars_A32, ...
    theta_o, psi, targets, params);
R_hi = vars_hi.Ry;

if R_lo * R_hi > 0
    error('Bisection fails: residual has same sign at bounds.');
end

tol = 1e-8;
maxit = 200;

for it = 1:maxit

    yl_mid = 0.5 * (yl_lo + yl_hi);
    vars_mid  = yl_residual(yl_mid, vars_A31, vars_A32, ...
    theta_o, psi, targets, params);
    R_mid = vars_mid.Ry;
    
    if abs(R_mid) < tol || abs(yl_hi - yl_lo) < tol
        vars = vars_mid;
        vars.yl = yl_mid;
        return
    end

    if R_lo * R_mid < 0
        yl_hi = yl_mid;
        R_hi  = R_mid;
    else
        yl_lo = yl_mid;
        R_lo  = R_mid;
    end


end % end yl loop 

error('Bisection for yl did not converge');

end


