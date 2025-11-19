/*
  Galí (2015) – Tabla 4.1 (Taylor rule, CURRENT, linear).
  Compatible con Dynare-Julia. Los valores por defecto se pueden
  sobreescribir con -D (sin espacios):
    -DSHOCKCASE=1|2   (1: Tecnología, 2: Demanda)
    -DPHI_PI=1.5
    -DPHI_Y=0.125
*/

// ===== Macros con fallback (permiten -D) =====
@#if !defined(SHOCKCASE)
  @#define SHOCKCASE 1
@#endif
@#if !defined(PHI_PI)
  @#define PHI_PI 1.5
@#endif
@#if !defined(PHI_Y)
  @#define PHI_Y 0.125
@#endif

// ===== Variables =====
var pi y_gap y_nat y yhat i r_nat a z;
varexo eps_a eps_z;

// ===== Parámetros =====
parameters betta siggma varphi alppha epsilon theta rho_a rho_z phi_pi phi_y;
parameters Omega psi_n_ya lambda kappa;

// --- Calibración base (Galí)
betta   = 0.99;
siggma  = 1;
varphi  = 5;
alppha  = 0.25;
epsilon = 9;
theta   = 0.75;

rho_a   = 0.90;
rho_z   = 0.50;

// --- Coeficientes de la regla (desde macros)
phi_pi  = @{PHI_PI};
phi_y   = @{PHI_Y};

// --- Parámetros compuestos
Omega     = (1-alppha)/(1-alppha+alppha*epsilon);
psi_n_ya  = (1+varphi)/(siggma*(1-alppha)+varphi+alppha);
lambda    = ((1-theta)*(1-betta*theta)/theta)*Omega;
kappa     = lambda*(siggma + (varphi+alppha)/(1-alppha));

// ===== Modelo lineal =====
model(linear);
  // NK Phillips curve
  pi = betta*pi(+1) + kappa*y_gap;

  // IS en brecha (gap)
  y_gap = -(1/siggma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Tasa natural y producto natural
  r_nat = -siggma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya*a;

  // Producto total y desvío
  y    = y_nat + y_gap;
  yhat = y;                        // en lineal, y ya es desvío

  // Regla de Taylor (corriente) con output gap
  i = phi_pi*pi + phi_y*y_gap;

  // Procesos de shocks
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;
end;

steady;
check;

// ===== Bloque de shocks (1: Tec.; 2: Demanda) =====
@#if SHOCKCASE==1
  shocks;
    var eps_a; stderr 1;
    var eps_z; stderr 0;
  end;
@#elseif SHOCKCASE==2
  shocks;
    var eps_a; stderr 0;
    var eps_z; stderr 1;
  end;
@#else
  shocks;
    var eps_a; stderr 1;
    var eps_z; stderr 0;
  end;
@#endif

// Momentos teóricos (Julia extrae sigmas de la salida)
stoch_simul(order=1, irf=0, nograph) pi y_gap y yhat i r_nat a z;
