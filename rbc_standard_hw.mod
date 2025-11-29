// Hansen & Wright (1992) — Standard RBC (Dynare/Dynare.jl)
// Solve at order=1. Big Monte Carlo + HP filter will be done in Julia.
//
// Variables: y output, c consumption, i investment, k capital, h hours, l leisure, z technology
// Shock: e ~ N(0,1), scaled by sigma_e inside z law of motion.

var c i y k h l z;
varexo e;

parameters beta delta theta rho sigma_e A;

// Calibration (as in your write-up)
beta    = 0.99;
delta   = 0.025;
theta   = 0.36;
rho     = 0.95;
sigma_e = 0.007;

// Chosen to target h_ss = 1/3 (your computation)
A       = 1.7214;

model;
    // identities/technology
    l = 1 - h;
    y = exp(z)*k^theta*h^(1-theta);

    // feasibility + laws of motion
    c + i = y;
    k(+1) = (1-delta)*k + i;
    z     = rho*z(-1) + sigma_e*e;

    // Euler
    1/c = beta*(1/c(+1))*( theta*exp(z(+1))*k(+1)^(theta-1)*h(+1)^(1-theta) + 1 - delta );

    // Intra-temporal (consumption vs leisure)
    A/l = (1/c)*(1-theta)*exp(z)*k^theta*h^(-theta);
end;

// Deterministic steady state (z=0, h=1/3 target)
steady_state_model;
    z = 0;
    h = 1/3;
    l = 1 - h;

    // From steady-state Euler:
    // 1 = beta*( theta*(k/h)^(theta-1) + 1 - delta )
    kh = ((1/beta - 1 + delta)/theta)^(1/(theta-1));   // k/h
    k  = kh*h;

    y = exp(z)*k^theta*h^(1-theta);
    i = delta*k;
    c = y - i;
end;

steady;
check;
resid(1);

// e has unit stdev; sigma_e is already inside the z equation
shocks;
    var e; stderr 1;
end;

// Solve only (no Dynare simulation here)
stoch_simul(order=1, irf=0, nograph, nodisplay);
