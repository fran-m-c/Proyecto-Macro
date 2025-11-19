/*
  Galí (2015) – Tabla 4.1 (Regla de Taylor, CURRENT).
  Versión “Julia-safe”: sin scripting MATLAB; la tabla se arma en Julia.

  Defines por línea de comando (sin espacios):
    -DSHOCKCASE=1   -> Technology shock
    -DSHOCKCASE=2   -> Demand shock
    -DPHI_PI=1.5
    -DPHI_Y=0.125
*/

@#define SHOCKCASE=1     // 1=TECH, 2=DEMAND
@#define PHI_PI=1.5
@#define PHI_Y=0.125

var pi y_gap y_nat y yhat i r_nat a z;
varexo eps_a eps_z;

parameters betta siggma varphi alppha epsilon theta rho_a rho_z phi_pi phi_y;
parameters Omega psi_n_ya lambda kappa;

// --- Calibración (Galí 2015)
betta  = 0.99;
siggma = 1;
varphi = 5;
alppha = 1/4;
epsilon= 9;
theta  = 3/4;

rho_a  = 0.90;
rho_z  = 0.50;

// --- Coeficientes de la regla (desde macros)
phi_pi = @{PHI_PI};
phi_y  = @{PHI_Y};

// --- Parámetros compuestos
Omega     = (1-alppha)/(1-alppha+alppha*epsilon);
psi_n_ya  = (1+varphi)/(siggma*(1-alppha)+varphi+alppha);
lambda    = (1-theta)*(1-betta*theta)/theta * Omega;
kappa     = lambda*(siggma + (varphi+alppha)/(1-alppha));

// --- Modelo lineal (desviaciones)
model(linear);
  // NK Phillips curve
  pi = betta*pi(+1) + kappa*y_gap;

  // IS en gaps
  y_gap = -(1/siggma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Natural rate y output natural
  r_nat = -siggma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya * a;

  // Output total, desvío y regla de Taylor (corriente)
  y    = y_nat + y_gap;
  yhat = y - steady_state(y);   // en lineal: yhat = y
  i    = phi_pi*pi + phi_y*y_gap;   // Taylor "current" sobre el GAP

  // Choques
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;
end;

steady; check;

// --- Bloque de choques (1=TECH, 2=DEMAND)
@#if SHOCKCASE==1
  shocks; var eps_a = 1; var eps_z = 0; end;
@#elseif SHOCKCASE==2
  shocks; var eps_a = 0; var eps_z = 1; end;
@#else
  shocks; var eps_a = 1; var eps_z = 0; end;
@#endif

// Momentos teóricos (leer desde Julia)
stoch_simul(order=1, irf=0, nograph);
