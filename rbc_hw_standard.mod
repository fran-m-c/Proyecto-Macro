var c i y k h z;
varexo eps;

parameters beta delta theta rho sigma A;

beta  = 0.99;
delta = 0.025;
theta = 0.36;
rho   = 0.95;
sigma = 0.007;

A = 1.7214;

model;

y = exp(z)*k^theta*h^(1-theta);
y = c + i;

k = (1-delta)*k(-1) + i;

z = rho*z(-1) + eps;

1/c = beta * (1/c(+1))*(exp(z(+1))*theta*k(+1)^(theta-1)*h(+1)^(1-theta) + 1-delta);

A/(1-h) = (1/c)*exp(z)*(1-theta)*k^theta*h^(-theta);

end;

shocks;
    var eps; stderr sigma;
end;

steady;
check;

stoch_simul(order=1, irf=0, nograph);

save_params_and_steady_state;
write_latex_dynamic_model;
