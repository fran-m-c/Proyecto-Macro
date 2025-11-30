var c i y k h z;
varexo e_z;

parameters beta delta theta rho A;

beta  = 0.99;
delta = 0.025;
theta = 0.36;
rho   = 0.95;
A     = 1.7214;

model;
    y = exp(z)*k(-1)^theta*h^(1-theta);
    c + i = y;
    k = (1-delta)*k(-1) + i;
    z = rho*z(-1) + e_z;

    1/c = beta*(1/c(+1))*(theta*exp(z(+1))*k^(theta-1)*h(+1)^(1-theta) + 1 - delta);
    A/(1-h) = (1/c) * (1-theta) * y / h;
end;

initval;
    z = 0;
    h = 1/3;
    k = 12.6629;
    y = 1.2347;
    i = 0.3166;
    c = 0.918096;
end;

shocks;
    var e_z; stderr 1;
end;

steady;
check;

stoch_simul(order=1, irf=0, noprint, nograph);
