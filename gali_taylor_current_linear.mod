/*
 * Gali (2015) – Table 4.1: Evaluation of Simple Rules (Taylor Rule)
 * LINEAR model in deviations from steady state.
 * This file computes both panels: Technology shocks and Demand shocks,
 * over the four (phi_pi, phi_y) configurations used in the book:
 *   (1.5, 0.125), (1.5, 0), (5, 0), (1, 1.5)  <-- note: last column is phi_y=1 (output), phi_pi=1.5
 * Outputs a dyntable with sigma(y), sigma(tilde y), sigma(pi), and welfare L.
 */

var pi          // inflation
    y_gap       // output gap \tilde y
    y_nat       // natural output
    y           // output
    yhat        // y - steady_state(y)
    i           // nominal interest rate
    r_nat       // natural real rate
    a           // technology process
    z;          // preference (demand) process

varexo eps_a eps_z;

parameters alppha betta rho_a rho_z siggma varphi epsilon theta phi_pi phi_y;
parameters Omega psi_n_ya lambda kappa; // composite parameters

// ---------------- Calibration (quarterly, as in Gali 2015)
siggma = 1;
varphi = 5;
alppha = 1/4;
epsilon= 9;
theta  = 3/4;
betta  = 0.99;

rho_a  = 0.90;
rho_z  = 0.50;

// Default values (will be overwritten in loops)
phi_pi = 1.5;
phi_y  = 0.125;

// ---------------- Composite parameters (Gali pp. 60–63)
Omega     = (1-alppha)/(1-alppha+alppha*epsilon);
psi_n_ya  = (1+varphi)/(siggma*(1-alppha)+varphi+alppha);
lambda    = (1-theta)*(1-betta*theta)/theta * Omega;
kappa     = lambda*(siggma + (varphi+alppha)/(1-alppha));

// ---------------- Linear model in deviations
model(linear);
  // NKPC
  pi = betta*pi(+1) + kappa*y_gap;

  // Dynamic IS (gap form)
  y_gap = -(1/siggma) * ( i - pi(+1) - r_nat ) + y_gap(+1);

  // Natural rate and natural output
  r_nat = -siggma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya*a;

  // Output identity and yhat
  y = y_nat + y_gap;
  yhat = y - steady_state(y); // = y in linear deviations, but keep it explicit

  // Taylor rule (CURRENT values, as in the book’s simple rule)
  i = phi_pi*pi + phi_y*yhat;

  // Shock processes
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;
end;

// ---------------- Solve base (get var list)
steady; check;
// run once (no printing) only to initialize internal structures
stoch_simul(order=1, irf=0, noprint) y y_gap pi;

// ---------------- Variable positions
y_pos     = strmatch("y",     var_list_, "exact");
y_gap_pos = strmatch("y_gap", var_list_, "exact");
pi_pos    = strmatch("pi",    var_list_, "exact");

// ---------------- Helper for welfare L
// L = 1/2 * [ (sigma + (varphi+alpha)/(1-alpha)) Var(y_gap) + (epsilon/lambda) Var(pi) ] / 100
// (the /100 scaling matches the convention used in the Dynare reference implementation)
W_weight = siggma + (varphi+alppha)/(1-alppha);

// ---------------- PANEL A: Technology shocks (Table 4.1 left block)
phi_pi_vec = [1.5, 1.5, 5.0, 1.5];
phi_y_vec  = [0.125, 0.0, 0.0, 1.0];

variance_y_gap = NaN(1, length(phi_pi_vec));
variance_y     = NaN(1, length(phi_pi_vec));
variance_pi    = NaN(1, length(phi_pi_vec));
Lval           = NaN(1, length(phi_pi_vec));

// set shocks: Technology only
shocks; var eps_a = 1^2; var eps_z = 0; end;

for ii = 1:length(phi_pi_vec)
  set_param_value("phi_pi", phi_pi_vec(ii));
  set_param_value("phi_y",  phi_y_vec(ii));

  [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); // re-solve

  if ~info(1)
    // read variances
    variance_y_gap(ii) = oo_.var(y_gap_pos, y_gap_pos);
    variance_y(ii)     = oo_.var(y_pos,     y_pos);
    variance_pi(ii)    = oo_.var(pi_pos,    pi_pos);
    // welfare
    Lval(ii) = 0.5 * ( W_weight*variance_y_gap(ii) + epsilon/lambda*variance_pi(ii) ) / 100;
  end
end

labels  = {"phi_pi"; "phi_y"; "sigma(y)"; "sigma(tilde y)"; "sigma(pi)"; "L"};
headers = {" ", " ", " ", " "};
values  = [phi_pi_vec; phi_y_vec; sqrt(variance_y); sqrt(variance_y_gap); sqrt(variance_pi); Lval];

// print panel A
options_.noprint = 0;
dyntable(options_, "Technology", headers, labels, values, size(labels,2)+2, 4, 3);
options_.noprint = 1;

// ---------------- PANEL B: Demand shocks (Table 4.1 right block)
variance_y_gap = NaN(1, length(phi_pi_vec));
variance_y     = NaN(1, length(phi_pi_vec));
variance_pi    = NaN(1, length(phi_pi_vec));
Lval           = NaN(1, length(phi_pi_vec));

// set shocks: Demand only
shocks; var eps_a = 0; var eps_z = 1^2; end;

for ii = 1:length(phi_pi_vec)
  set_param_value("phi_pi", phi_pi_vec(ii));
  set_param_value("phi_y",  phi_y_vec(ii));

  [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);

  if ~info(1)
    variance_y_gap(ii) = oo_.var(y_gap_pos, y_gap_pos);
    variance_y(ii)     = oo_.var(y_pos,     y_pos);
    variance_pi(ii)    = oo_.var(pi_pos,    pi_pos);
    Lval(ii) = 0.5 * ( W_weight*variance_y_gap(ii) + epsilon/lambda*variance_pi(ii) ) / 100;
  end
end

values = [phi_pi_vec; phi_y_vec; sqrt(variance_y); sqrt(variance_y_gap); sqrt(variance_pi); Lval];

// print panel B
options_.noprint = 0;
dyntable(options_, "Demand", headers, labels, values, size(labels,2)+2, 4, 3);
options_.noprint = 1;
