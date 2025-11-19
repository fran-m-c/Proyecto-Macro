/*
  NK lineal (Gali 2015) + Regla de Taylor "corriente":
    i_t = phi_pi * pi_t + phi_y * y~_t
  Conmutador de choques: TECH / DEMAND / BOTH
*/

@#define SHOCKCASE = "BOTH"   // <-- cambia a "TECH" o "DEMAND"

var pi y_gap y_nat y i r_nat a z;
varexo eps_a eps_z;

parameters beta sigma varphi alpha epsilon theta rho_a rho_z phi_pi phi_y;
parameters Omega psi_n_ya lambda kappa;

// --- Calibración (Gali 2015)
beta   = 0.99;
sigma  = 1;
varphi = 5;
alpha  = 0.25;
epsilon= 9;
theta  = 3/4;

rho_a  = 0.90;
rho_z  = 0.50;

phi_pi = 1.5;
phi_y  = 0.125;

// --- Comp. (p.60–63)
Omega     = (1-alpha)/(1-alpha+alpha*epsilon);
psi_n_ya  = (1+varphi)/(sigma*(1-alpha)+varphi+alpha);
lambda    = (1-theta)*(1-beta*theta)/theta * Omega;
kappa     = lambda*(sigma + (varphi+alpha)/(1-alpha));

model(linear);
  // NKPC
  pi = beta*pi(+1) + kappa*y_gap;

  // DIS
  y_gap = -(1/sigma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Natural rate / natural output
  r_nat = -sigma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya * a;

  // Output total
  y = y_nat + y_gap;

  // Regla de Taylor (corriente)
  i = phi_pi*pi + phi_y*y_gap;

  // Choques
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;
end;

steady; check;

// --- Selección de choques para momentos/IRFs
@#if SHOCKCASE == "TECH"
  shocks;
    var eps_a = 1;
    var eps_z = 0;
  end;
@#elseif SHOCKCASE == "DEMAND"
  shocks;
    var eps_a = 0;
    var eps_z = 1;
  end;
@#else
  shocks;
    var eps_a = 1;
    var eps_z = 1;
  end;
@#endif

// Teóricos + IRFs (ajusta irf si quieres)
stoch_simul(order=1, irf=15, nograph);
