function vars = a32_ownership_other(yh, xh, psi, theta_o, targets, params)

%========================================================================%
% Appendix A.3.2: Other ownership-market variables
%
% Inputs:
%   yh, xh      = thresholds from Appendix A.3.1
%   psi,theta_o = guesses for equilibrium vars
%   targets     = Table 3 calibration targets
%   params      = Table 4 structural parameters
%
% Outputs:
% vars.pi_h   = pi_h;         % Probability buyer finds acceptable match
% vars.Sigma_h = Sigma_h;     % Surplus from matching in ownership market
% vars.s_o     = s_o;         % Seller meeting rate
% vars.m_h     = m_h;         % Moving hazard for homeowners
% vars.q_h     = q_h;         % Owner-occupier stock
% vars.b_h     = b_h;         % Buyer stock
% vars.u_o     = u_o;         % Seller stock
% vars.vo      = vo;          % Buyer meeting rate
% vars.G       = G;           % Government Revenue

%========================================================================%

%% Unpack parameters

zeta_h   = params.zeta_h;
lambda_h = params.lambda_h;
delta_h  = params.delta_h;

alpha_h  = params.alpha_h;
r        = params.r;
rho      = params.rho;

omega_h  = params.omega_h;
omega_i  = params.omega_i;

nu_o     = params.nu_o;
eta_o    = params.eta_o;

tau_h    = targets.tau_h;


% Market tightness determines meeting rate
vo = nu_o * theta_o^(-eta_o);

%% ----------------------------------------------------------
% 1. Transaction probability π_h = (z_h / y_h)^{λ_h}
%% ----------------------------------------------------------

pi_h = (zeta_h / yh)^lambda_h;

%% ----------------------------------------------------------
% 2. Surplus term Σ_h(yh,xh) = (A.51)
%% ----------------------------------------------------------

coeff = zeta_h^lambda_h / ((1 + tau_h*omega_h) * (r + rho + alpha_h) * (lambda_h - 1));

term1 = yh^(1 - lambda_h);
term2 = alpha_h * delta_h^lambda_h * xh^(1 - lambda_h) / ...
       (r + rho + alpha_h*(1 - delta_h^lambda_h));

Sigma_h = coeff * vo * (term1 + term2);

%% ----------------------------------------------------------
% 3. Seller-side meeting rate s_o = θ_o v_o [ ψ + (1-ψ)π_h ]   (A.56)
%% ----------------------------------------------------------

s_o = theta_o * vo * (psi + (1 - psi)*pi_h);

%% ----------------------------------------------------------
% 4. Buyer moving hazard m_h   (A.78)
%% ----------------------------------------------------------

m_h = alpha_h * (1 - (delta_h^(lambda_h) * (xh^(lambda_h - 1)) / (yh^(lambda_h - 1))));

%% ----------------------------------------------------------
% 5. Stocks from steady-state conditions (A.92)
%% ----------------------------------------------------------

Tmh = targets.Tmh;
Tbh = targets.Tbh;
Tso = targets.Tso;

n    = targets.n;
h    = targets.h;
i_sh = targets.i_share;

% Owner-occupiers q_h
q_h = n * h * Tmh / (Tmh + Tbh);

% Buyers b_h
b_h = (Tbh / Tmh) * q_h;

% Sellers u_o
u_o = Tso * q_h / ((1 - i_sh) * (Tmh + Tbh));

% ----------------------------------------------------------
% 6. Government revenue G from land transfer taxes (A.108)
% ----------------------------------------------------------

tau = targets.tau_i;     % same as tau_h in baseline calibration
G   = tau * s_o * u_o;

%% Pack results
vars = struct();
vars.pi_h   = pi_h;         % Probability buyer finds acceptable match
vars.Sigma_h = Sigma_h;     % Surplus from matching in ownership market
vars.s_o     = s_o;         % Seller meeting rate
vars.m_h     = m_h;         % Moving hazard for homeowners
vars.q_h     = q_h;         % Owner-occupier stock
vars.b_h     = b_h;         % Buyer stock
vars.u_o     = u_o;         % Seller stock
vars.vo      = vo;          % Buyer meeting rate
vars.G       = G;           % Government Revenue

end