// rbc_standard_hw_simple.mod
// Standard RBC (Hansen & Wright-style calibration)

var c i y k h l z;
varexo eps;

parameters beta delta theta rho sigma_e A;

beta    = 0.99;
delta   = 0.025;
theta   = 0.36;
rho     = 0.95;
sigma_e = 0.007;

// A se calibra endógenamente en steady_state_model para que h = 1/3

model;
    // Producción
    y = exp(z)*k^theta*h^(1-theta);

    // Restricción de recursos
    y = c + i;

    // Ocio
    l = 1 - h;

    // Acumulación de capital
    k(+1) = (1-delta)*k + i;

    // Euler
    1/c = beta*(1/c(+1))*(exp(z(+1))*theta*k(+1)^(theta-1)*h(+1)^(1-theta) + 1 - delta);

    // Intratemporal consumo-ocio
    A/(1-h) = (1/c)*exp(z)*(1-theta)*k^theta*h^(-theta);

    // Proceso tecnológico
    z = rho*z(-1) + eps;
end;

steady_state_model;
    z = 0;
    h = 1/3;
    l = 1 - h;

    k = ((1/beta - 1 + delta)/(theta*h^(1-theta)))^(1/(theta-1));
    y = exp(z)*k^theta*h^(1-theta);
    i = delta*k;
    c = y - i;

    // Calibra A para que h = 1/3 sea consistente con la FOC intratemporal
    A = ((1-theta)*exp(z)*k^theta*h^(-theta))*(1-h)/c;
end;

shocks;
    var eps; stderr sigma_e;
end;

steady;
check;

// Solo para que Dynare compute la solución (no usamos sus simulaciones)
stoch_simul(order=1, irf=0, periods=200, nograph);
