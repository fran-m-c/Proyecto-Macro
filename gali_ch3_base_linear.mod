/*
  Gali (2015), Cap. 3 — Modelo NK básico (lineal) en desviaciones
  Cierre "base" solo para funciones de política: i_t = 0.
  Para reglas de Taylor usar los otros .mod.
*/

var pi y_gap y_nat y i r_nat a z;   // inflación, brecha, output natural/total, tasa nominal, r^nat, choques
varexo eps_a eps_z;                 // choques de tecnología (a) y preferencias (z)

parameters beta sigma varphi alpha epsilon theta rho_a rho_z;
parameters Omega psi_n_ya lambda kappa; // compuestos

// ---- Calibración (trimestral, Galí 2015)
beta   = 0.99;
sigma  = 1;
varphi = 5;
alpha  = 0.25;
epsilon= 9;
theta  = 3/4;

rho_a  = 0.90;
rho_z  = 0.50;

// ---- Parámetros compuestos (Galí, pp. 60–63)
Omega     = (1-alpha)/(1-alpha+alpha*epsilon);
psi_n_ya  = (1+varphi)/(sigma*(1-alpha)+varphi+alpha);
lambda    = (1-theta)*(1-beta*theta)/theta * Omega;
kappa     = lambda*(sigma + (varphi+alpha)/(1-alpha));

// ---- Modelo (lineal)
model(linear);
  // NKPC
  pi = beta*pi(+1) + kappa*y_gap;

  // DIS (gap)
  y_gap = -(1/sigma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Tasa natural y output natural
  r_nat = -sigma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya * a;

  // Output total y definición de gap
  y = y_nat + y_gap;

  // Procesos AR(1)
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;

  // Cierre base: tasa nominal fija en su SS (desviación = 0)
  i = 0;
end;

// ---- Estacionario (todo 0 por linealidad) y chequeos
steady; check;

// ---- Choques (ajusta para casos "solo tecnología" o "solo demanda")
shocks;
  var eps_a = 1;     // varianza unidad
  var eps_z = 1;
end;

// ---- Política de simulación: imprime funciones de decisión
stoch_simul(order=1, irf=0, nograph);
