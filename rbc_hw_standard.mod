// Hansen & Wright (1992) – Modelo RBC estándar, versión mínima
// Archivo: rbc_hw_standard.mod

var c i y k h z;
varexo eps;

parameters beta delta theta rho sigma_eps A;

// ---- Parámetros (tus valores) ----
beta      = 0.99;
delta     = 0.025;
theta     = 0.36;
rho       = 0.95;
sigma_eps = 0.007;
A         = 1.7214; // calculado a partir de h=1/3

// ---- Modelo ----
model;

// Producción
y = exp(z)*k^theta*h^(1-theta);

// Restricción de recursos
y = c + i;

// Acumulación de capital
k(+1) = (1-delta)*k + i;

// Proceso de productividad
z = rho*z(-1) + eps;

// Euler
1/c = beta*(1/c(+1))*(exp(z(+1))*theta*k(+1)^(theta-1)*h(+1)^(1-theta) + 1-delta);

// FOC ocio-trabajo (intratemporal)
A/(1-h) = (1/c)*exp(z)*(1-theta)*k^theta*h^(-theta);

end;

// ---- Shock ----
shocks;
var eps; stderr sigma_eps;
end;

// ---- Steady state numérico sencillo (solo punto de partida) ----
initval;
z = 0;

// h=1/3
h = 1/3;

// k/h ratio desde Euler en estado estacionario:
// 1 = beta*(theta*(k/h)^(theta-1) + 1-delta)
k = 12.6631;   // tu valor
y = exp(z)*k^theta*h^(1-theta);
i = delta*k;
c = y - i;
end;

steady;

// No necesito que Dynare simule nada aquí: lo hará Julia.
// Pero para obligar a Dynare a calcular la solución lineal:
stoch_simul(order=1, irf=0, nograph);

