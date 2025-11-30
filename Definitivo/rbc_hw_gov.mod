var c i y k h z g;
varexo e_z e_g;

parameters beta delta theta rho A lamg ghat;

beta  = 0.99;
delta = 0.025;
theta = 0.36;
rho   = 0.95;

lamg  = 0.96;
ghat  = 0.271634;     // = 0.22 * y_ss (con y_ss≈1.2347)
A     = 2.4445;       // calibrado con h=1/3 y g/y=0.22

model;
    y = exp(z)*k(-1)^theta*h^(1-theta);
    c + i + g = y;

    k = (1-delta)*k(-1) + i;
    z = rho*z(-1) + e_z;

    // AR(1) en logs para g_t (en niveles en el modelo)
    log(g) = (1-lamg)*log(ghat) + lamg*log(g(-1)) + e_g;

    1/c = beta*(1/c(+1))*(theta*exp(z(+1))*k^(theta-1)*h(+1)^(1-theta) + 1 - delta);
    A/(1-h) = (1/c) * (1-theta) * y / h;
end;

initval;
    z = 0;
    h = 1/3;
    k = 12.6629;
    y = 1.2347;
    i = 0.3166;
    g = 0.271634;
    c = 0.6465;
end;

shocks;
    var e_z; stderr 1;
    var e_g; stderr 1;
end;

steady;
check;

stoch_simul(order=1, irf=0, noprint, nograph);
