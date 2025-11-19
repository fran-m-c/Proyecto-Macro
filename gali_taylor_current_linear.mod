/*
  Gali (2015) Table 4.1 – Taylor rule (CURRENT):
    i_t = phi_pi * pi_t + phi_y * y_gap_t
  Archivo “Julia-safe”: sin scripting MATLAB; todo post-proceso se hace en Julia.

  Usa defines desde Julia:
    defines = Dict("SHOCKCASE"=>"TECH" o "DEMAND",
                   "PHI_PI"=>"1.5", "PHI_Y"=>"0.125")
*/

@#define SHOCKCASE = "TECH"   // <- se sobreescribe desde Julia
@#define PHI_PI    = 1.5
@#define PHI_Y     = 0.125

var pi y_gap y_nat y i r_nat a z;
varexo eps_a eps_z;

parameters betta siggma varphi alppha epsilon theta rho_a rho_z phi_pi phi_y;
parameters Omega psi_n_ya lambda kappa;

// -------- Calibración (Gali 2015)
betta  = 0.99;
siggma = 1;
varphi = 5;
alppha = 1/4;
epsilon= 9;
theta  = 3/4;

rho_a  = 0.90;
rho_z  = 0.50;

// (posibles overrides desde Julia)
phi_pi = @{PHI_PI};
phi_y  = @{PHI_Y};

// -------- Parámetros compuestos
Omega     = (1-alppha)/(1-alppha+alppha*epsilon);
psi_n_ya  = (1+varphi)/(siggma*(1-alppha)+varphi+alppha);
lambda    = (1-theta)*(1-betta*theta)/theta * Omega;
kappa     = lambda*(siggma + (varphi+alppha)/(1-alppha));

// -------- Modelo lineal
model(linear);
  // NKPC
  pi = betta*pi(+1) + kappa*y_gap;

  // IS (gap)
  y_gap = -(1/siggma)*( i - pi(+1) - r_nat ) + y_gap(+1);

  // Tasas/producción “naturales”
  r_nat = -siggma*psi_n_ya*(1-rho_a)*a + (1-rho_z)*z;
  y_nat = psi_n_ya * a;

  // Output total y regla de Taylor (corriente)
  y = y_nat + y_gap;
  i = phi_pi*pi + phi_y*y_gap;

  // Choques
  a = rho_a*a(-1) + eps_a;
  z = rho_z*z(-1) + eps_z;
end;

steady; check;

// -------- Escoge el bloque de choques
@#if SHOCKCASE == "TECH"
  shocks; var eps_a = 1; var eps_z = 0; end;
@#elseif SHOCKCASE == "DEMAND"
  shocks; var eps_a = 0; var eps_z = 1; end;
@#else
  @#error Bad SHOCKCASE (use TECH or DEMAND)
@#endif

// Solo momentos teóricos (sin IRFs)
stoch_simul(order=1, irf=0, nograph);
