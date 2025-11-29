var c i y k h z;
varexo eps;

parameters beta delta theta rho sigma A;

beta  = 0.99;
delta = 0.025;
theta = 0.36;
rho   = 0.95;
sigma = 0.007;

// A calibrado externamente según tu derivación
A     = 1.7214;

model;

// Tecnología
y = exp(z)*k^theta*h^(1-theta);

// Recursos
y = c + i;

// Acumulación de capital
k = (1-delta)*k(-1) + i;

// Shock
z = rho*z(-1) + eps;

// Euler
1/c = beta * (1/c(+1))*(exp(z(+1))*theta*k(+1)^(theta-1)*h(+1)^(1-theta) + 1-delta);

// Intratemporal
A/(1-h) = (1/c)*exp(z)*(1-theta)*k^theta*h^(-theta);

end;

shocks;
    var eps; stderr sigma;
end;

steady;
check;

// Esto obliga a Dynare a calcular la solución lineal
stoch_simul(order=1, irf=0, nograph);
