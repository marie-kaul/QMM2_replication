function vars = a32_ownership_other(vars_A31, psi, theta_o, targets, params)

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
% vars.s_o     = s_o;         % Seller meeting rate
% vars.m_h     = m_h;         % Moving hazard for homeowners
% vars.q_h     = q_h;         % Owner-occupier stock
% vars.b_h     = b_h;         % Buyer stock
% vars.u_o     = u_o;         % Seller stock
% vars.vo      = vo;          % Buyer meeting rate
% vars.G       = G;           % Government Revenue

%========================================================================%

%% Unpack parameters and variables

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

yh       = vars_A31.yh;
xh       = vars_A31.xh;
Sigma_h  = vars_A31.Sigma_h;


% Market tightness determines meeting rate
vo = nu_o * theta_o^(-eta_o);

%% ----------------------------------------------------------
% 1. Transaction probability π_h = (z_h / y_h)^{λ_h}
%% ----------------------------------------------------------

pi_h = (zeta_h / yh)^lambda_h;

%% ----------------------------------------------------------
% 2. Seller-side meeting rate s_o = θ_o v_o [ ψ + (1-ψ)π_h ]   (A.56)
%% ----------------------------------------------------------

s_o = theta_o * vo * (psi + (1 - psi)*pi_h);

%% ----------------------------------------------------------
% 3. Buyer moving hazard m_h   (A.58)
%% ----------------------------------------------------------

delta_lam = delta_h^lambda_h;
ratio_lam = (yh / xh)^lambda_h;

num = rho + alpha_h * (1 - delta_lam) ...
      - rho * delta_lam * ratio_lam;

den = rho + alpha_h * (1 - delta_lam) ...
      + alpha_h * delta_lam * ratio_lam;

m_h = alpha_h * (num / den);

%% ----------------------------------------------------------
% 4. Stocks from steady-state conditions (A.56)–(A.59)
%% ----------------------------------------------------------

% u_o from (A.59)
u_o = 1 / ( ...
      1 ...
    + (1 - psi) * s_o / (m_h + rho) ...
    + psi * s_o / params.rho_l );

% q_h from (A.59)
q_h = ((1 - psi) * s_o / (m_h + rho)) * u_o;

% buyers
b_h = (1 - psi) * theta_o * u_o;
b_i = psi * theta_o * u_o;


%% ----------------------------------------------------------
% 5. Government revenue G from land transfer taxes 
%% ----------------------------------------------------------

% Transaction shares
i = psi / (psi + (1-psi)*pi_h);

% Total sales flow
S_total = s_o * u_o;
S_h = (1 - i) * S_total;
S_i = i * S_total;

% Prices
Sigma_i = vars_A31.Sigma_i; 

Co = params.Co; D = params.D;

% A.54
Ph = ((r + theta_o*vo*(1-psi)*pi_h)/r) * (omega_h*Sigma_h/pi_h) ...
     + (theta_o*vo*psi*omega_i*Sigma_i)/r ...
     + Co - D/r;

% A.55
Pi = Co ...
     + (theta_o*vo*((1-psi)*omega_h*Sigma_h + psi*omega_i*Sigma_i) - D)/r ...
     + omega_i*Sigma_i;


G = targets.tau_h * Ph * S_h + targets.tau_i * Pi * S_i;

%% ----------------------------------------------------------
% 6. Pack Results
%% ----------------------------------------------------------

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
vars.b_i     = b_i;


end