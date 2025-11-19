/*
  Gali (2015) – Table 4.1 (Taylor rule, CURRENT).
  Versión “Julia-safe”: sin scripting MATLAB; la tabla se arma en Julia.

  Defines por línea de comando (sin espacios):
    -DSHOCKCASE=1   -> Technology shock
    -DSHOCKCASE=2   -> Demand shock
    -DPHI_PI=1.5
    -DPHI_Y=0.125
*/

/* ===== Macros con valores por defecto (sobrescribibles desde Julia) ===== */
@#ifndef SHOCKCASE
  @#define SHOCKCASE = 1     // 1 = Technology, 2 = Demand
@#endif

@#ifndef PHI_PI
  @#define PHI_PI = 1.5
@#endif

@#ifndef PHI_Y
  @#define PHI_Y = 0.125
@#endif

/* ===== Variables ===== */
var pi y_gap y_nat y yhat i r_nat a z;
varexo eps_a eps_z;

parameters betta siggma varphi alppha epsilon theta rho_a rho_z phi_pi phi_y;
parameters Omega psi_n_ya lambda kappa;

/* ===== Calibración (Galí 2015) ===== */
betta  = 0.99;
siggma = 1;
varphi = 5;
alppha = 1/4;
epsilon = 9;
theta   = 3/4;

rho_a = 0.90;
rho_z = 0.50;

/* ===== Coeficientes de la regla de Taylor (desde macros) ===== */
phi_pi = @{PHI_PI};
phi_y  = @{PHI_Y};

/* ===== Parámetros compuestos ===== */
Omega    = (1-alppha)/(1-alppha+alppha*epsilon);
psi_n_ya = (1+varphi)/(siggma*(1-alppha)+varphi+alppha);
lambda   = (1-theta)*(1-betta*theta)/theta * Omega;
kappa    = lambda*(siggma + (varphi+alppha)/(1-alppha));

/* ===== Modelo lineal (desviaciones respecto al estado estacionario) ===== */
model(linear);

  // Curva de Phillips NK
  pi = betta*pi(+1) + kappa*y_gap;

  // IS en brechas (output gap)
  y_gap = -(1/siggma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Tasa natural de interés y producto natural
  r_nat = -siggma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya * a;

  // Output total, desvío de output y regla de Taylor (corriente)
  y    = y_nat + y_gap;
  yhat = y - steady_state(y);    // en modelo lineal: yhat = y
  i    = phi_pi*pi + phi_y*yhat;

  // Procesos de choques
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;

end;

steady;
check;

/* ===== Bloque de choques (1 = Technology, 2 = Demand) ===== */
@#if SHOCKCASE==1
  // Solo choque tecnológico
  shocks;
    var eps_a = 1;
    var eps_z = 0;
  end;
@#elseif SHOCKCASE==2
  // Solo choque de demanda
  shocks;
    var eps_a = 0;
    var eps_z = 1;
  end;
@#else
  // Por defecto: tecnología
  shocks;
    var eps_a = 1;
    var eps_z = 0;
  end;
@#endif

/* ===== Momentos teóricos (leídos desde Julia) ===== */
stoch_simul(order=1, irf=0, nograph);
