/*
  Gali (2015) cap. 3 — Modelo NK básico (lineal, en desviaciones)
  Archivo "base" sin regla de política (cierre: i = 0).
*/

var pi y_gap y_nat y i r_nat a z;
varexo eps_a eps_z;

parameters beta sigma varphi alpha epsilon theta rho_a rho_z;
parameters Omega psi_n_ya lambda kappa;

// --- Calibración (trimestral, Gali 2015)
beta   = 0.99;
sigma  = 1;
varphi = 5;
alpha  = 0.25;
epsilon= 9;
theta  = 3/4;

rho_a  = 0.90;
rho_z  = 0.50;

// --- Parámetros compuestos
Omega     = (1-alpha)/(1-alpha+alpha*epsilon);
psi_n_ya  = (1+varphi)/(sigma*(1-alpha)+varphi+alpha);
lambda    = (1-theta)*(1-beta*theta)/theta * Omega;
kappa     = lambda*(sigma + (varphi+alpha)/(1-alpha));

// --- Modelo (lineal)
model(linear);
  // NKPC
  pi = beta*pi(+1) + kappa*y_gap;

  // DIS
  y_gap = -(1/sigma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Natural rate y output natural
  r_nat = -sigma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya * a;

  // Output total y definición de gap
  y = y_nat + y_gap;

  // Procesos AR(1)
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;

  // Cierre "base" (no determinante, solo para compilar el bloque)
  i = 0;
end;

steady; check;

// choques (solo para moverse; en este base no simules)
shocks;
  var eps_a = 1;
  var eps_z = 1;
end;

stoch_simul(order=1, irf=0, nograph); // imprime funciones si estuviera bien cerrado
